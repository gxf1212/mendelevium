---
title: "Deploy PostgreSQL Database and MinIO Object Storage: Complete Server Setup Guide"
date: "2025-05-26"
tags: [postgresql, minio, deployment, ubuntu, docker, nginx, reverse-proxy, database, object-storage]
---
# Ubuntu 22.04 服务器部署 PostgreSQL 数据库、MinIO 对象存储以及一个通过 Nginx 反向代理访问的 Docker化 Django 后端应用完整教程

**目标：** 部署 PostgreSQL 数据库、MinIO 对象存储以及一个通过 Nginx 反向代理访问的 Docker化 Django 后端应用。

**服务器 IP 定义 (请在脚本和配置中替换为您真实的服务器 IP)**： `SERVER_IP="123.45.6.78"` (示例 IP，请务必修改)

## **第 1 步：系统初始化与基础依赖安装**

首先，更新您的服务器并安装一些必要的工具。

```bash
# 更新系统包列表并升级现有包
sudo apt update && sudo apt upgrade -y

# 安装基础工具：ca-certificates, curl, gnupg, lsb-release 用于添加 Docker源，nginx 用于反向代理，git 用于拉取代码
sudo apt install -y ca-certificates curl gnupg lsb-release nginx git
```

## **第 2 步：安装 Docker CE 和 Docker Compose**

我们将使用 Docker 来容器化 MinIO 和 Django 应用。

```bash
# 1. 添加 Docker 官方 GPG 密钥
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

# 2. 设置 Docker APT 软件源
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

# 3. 安装 Docker CE (社区版), CLI, Containerd, 和 Docker Compose 插件
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# 4. （可选）配置国内 Docker 镜像加速器，以提高拉取镜像的速度
# 请根据您选择的云服务商或镜像源替换下面的地址
sudo mkdir -p /etc/docker
sudo tee /etc/docker/daemon.json <<EOF
{
  "registry-mirrors": [
    "https://docker.mirrors.ustc.edu.cn",
    "https://hub-mirror.c.163.com",
    "https://registry.docker-cn.com"
  ]
}
EOF
sudo systemctl daemon-reload
sudo systemctl restart docker

# 5. 将当前用户添加到 docker 组，这样执行 docker 命令时无需 sudo（需要重新登录或执行 newgrp docker 生效）
sudo usermod -aG docker $USER
echo "请重新登录或执行 'newgrp docker' 使 docker 组权限生效。"

# 6. 验证 Docker 是否安装成功
docker --version
docker compose version
```

## **第 3 步：安装 PostgreSQL 并配置远程访问**

Django 应用将使用 PostgreSQL作为数据库。

```bash
# 1. 安装 PostgreSQL 和相关工具
sudo apt install -y postgresql postgresql-contrib

# 2. 修改 PostgreSQL 配置以允许远程连接
# 编辑 postgresql.conf 文件，将 listen_addresses 从 'localhost' 改为 '*'
# 注意：您的 PostgreSQL 版本可能不同，请相应调整路径（例如 /etc/postgresql/16/main/）
# 您可以通过 `pg_lsclusters` 查看版本和路径
sudo sed -i "s/#listen_addresses = 'localhost'/listen_addresses = '*'/g" /etc/postgresql/$(pg_lsclusters | awk 'NR==2 {print $1}')/main/postgresql.conf

# 3. 修改 pg_hba.conf 文件以允许来自任何 IP 地址的 md5 密码认证连接
# 同样，注意 PostgreSQL 版本路径
echo "host    all             all             0.0.0.0/0               md5" | sudo tee -a /etc/postgresql/$(pg_lsclusters | awk 'NR==2 {print $1}')/main/pg_hba.conf

# 4. 重启 PostgreSQL 服务使配置生效
sudo systemctl restart postgresql

# 5. 设置 PostgreSQL 的 postgres 用户密码（重要！）
# 将 'YourSecurePostgresPassword!' 替换为您自己的强密码
sudo -u postgres psql -c "ALTER USER postgres WITH PASSWORD 'YourSecurePostgresPassword!';"
```

> **可视化操作与安全组说明 (PostgreSQL)**
>
> 1. **云服务器安全组**:
>
>    - 在您的云服务提供商（如阿里云、腾讯云、AWS）的控制台中，找到您的服务器实例对应的安全组（或防火墙规则）。
>    - 添加入站规则，允许来自您需要访问数据库的 IP 地址（或者为了开发方便暂时允许 `0.0.0.0/0`，但生产环境不推荐）访问 PostgreSQL 的默认端口 `5432/TCP`。
>
> 2. **数据库客户端连接**:
>
>    - 您可以使用图形化数据库管理工具（如 pgAdmin, DBeaver, Navicat 等）从您的本地计算机连接到服务器上的 PostgreSQL。
>    - 连接信息：
>      - **主机/服务器地址**: `YOUR_SERVER_IP` (例如 `123.45.6.78`)
>      - **端口**: `5432`
>      - **数据库**: 默认可以是 `postgres`
>      - **用户名**: `postgres`
>      - **密码**: 您在上面第 5 步设置的 `YourSecurePostgresPassword!`
>
> 3. **创建专用数据库和用户 (推荐)**: 虽然您可以使用 `postgres` 超级用户，但更安全的做法是为您的 Django 应用创建一个专用的数据库和用户。登录后，在 SQL 工具中执行：
>
>    ```bash
>    CREATE DATABASE myproject_db;
>    CREATE USER myproject_user WITH PASSWORD 'MyProjectSecurePassword!';
>    GRANT ALL PRIVILEGES ON DATABASE myproject_db TO myproject_user;
>    ALTER ROLE myproject_user CREATEDB; -- 可选，如果需要用户创建数据库
>    ```
>
>    之后在 Django 的 `settings.py` 中使用这些新的凭据。

## **第 4 步：部署 MinIO 对象存储 (使用 Docker)**

MinIO 将用于存储 Django 应用的媒体文件和静态文件。

```bash
wget https://dl.min.io/client/mc/release/linux-amd64/mc
chmod +x mc
sudo mv mc /usr/local/bin/

# 1. 创建 MinIO 数据存储目录
sudo mkdir -p /minio/data
sudo chmod -R 777 /minio/data # 临时给予宽松权限，生产环境应更精细控制

# 2. 使用 Docker 启动 MinIO 容器
# 将 'YourMinioAdminUser' 和 'YourMinioAdminPassword!' 替换为您自己的凭据
# 确保密码足够复杂（至少8位，包含大小写、数字、特殊字符）
docker run -d \
  --name minio \
  -p 9000:9000 \
  -p 9001:9001 \
  -v /minio/data:/data \
  -e "MINIO_ROOT_USER=admin" \
  -e "MINIO_ROOT_PASSWORD=YourSecureMinioPassword!" \
  quay.io/minio/minio:latest \
  server /data --console-address ":9001"
  
# 2. 设置别名（注意用单引号包裹密码）
mc alias set myminio http://123.45.6.78:9000 admin 'YourSecureMinioPassword!'
# 2. 修改密码
mc admin user password myminio admin 'NewSecurePass123!'
```

> **可视化操作与安全组说明 (MinIO)**
>
> 1. **云服务器安全组**:
>    - 开放 MinIO API 端口: `9000/TCP`
>    - 开放 MinIO 控制台端口: `9001/TCP`
> 2. **访问 MinIO 控制台**:
>    - 在浏览器中打开 `http://YOUR_SERVER_IP:9001` (例如 `http://123.45.6.78:9001`)。
>    - 使用您在 `docker run` 命令中设置的 `MINIO_ROOT_USER` (例如 `admin`) 和 `MINIO_ROOT_PASSWORD` (例如 `YourSecureMinioPassword!`) 登录。
> 3. **创建存储桶 (Buckets)**:
>    - 登录 MinIO 控制台后，点击 "Buckets" -> "Create Bucket"。
>    - 创建两个存储桶：
>      - `media` (用于存储用户上传的文件)
>      - `static` (用于存储 Django 的静态文件)
>    - **重要**: 为这两个存储桶设置访问策略 (Access Policy)。对于公开访问的静态文件和媒体文件，您可能需要将策略设置为 `public` 或 `readonly` (根据需求)。点击存储桶旁边的 "Manage" -> "Access Policy"，选择 "Add Policy"，然后选择 `readonly` 或 `download` (对于 `public`，可以直接在创建时设置)。更精细的权限控制请参考 MinIO 文档。

## **第 5 步：部署 Django 后端应用 (使用 Docker)**

### **5.1 准备 Django 项目代码**

```bash
# 1. 克隆您的 Django 项目代码 (替换为您的仓库地址)
# git clone https://github.com/your-username/your-backend-repo.git
# cd your-backend-repo
# 假设您已将代码上传到服务器的某个目录，例如 /srv/django-app
# cd /srv/django-app

# 2. 配置 Django settings.py (或通过环境变量传递)
#    确保您的 settings.py 文件中数据库和 MinIO 配置正确。
#    数据库示例 (使用您在 PostgreSQL 步骤中创建的用户和数据库):
#    DATABASES = {
#        'default': {
#            'ENGINE': 'django.db.backends.postgresql',
#            'NAME': 'myproject_db',
#            'USER': 'myproject_user',
#            'PASSWORD': 'MyProjectSecurePassword!',
#            'HOST': 'YOUR_SERVER_IP', # Django 容器需要能访问到宿主机的 PostgreSQL
#            'PORT': '5432',
#        }
#    }
#
#    MinIO (django-storages) 示例:
#    DEFAULT_FILE_STORAGE = 'storages.backends.s3boto3.S3Boto3Storage'
#    STATICFILES_STORAGE = 'storages.backends.s3boto3.S3Boto3Storage'
#    AWS_ACCESS_KEY_ID = 'admin' # MinIO Root User
#    AWS_SECRET_ACCESS_KEY = 'YourSecureMinioPassword!' # MinIO Root Password
#    AWS_STORAGE_BUCKET_NAME = 'media' # 默认文件存储桶
#    AWS_S3_ENDPOINT_URL = 'http://YOUR_SERVER_IP:9000' # MinIO API 地址
#    AWS_S3_OBJECT_PARAMETERS = { 'CacheControl': 'max-age=86400', }
#    AWS_DEFAULT_ACL = None # 或 'public-read' 根据需求
#    AWS_S3_USE_SSL = False # 如果 MinIO 没有配置 SSL
#    AWS_S3_VERIFY = True # 如果 MinIO 使用自签名证书，可能需要设为 False 或提供证书路径
#    AWS_S3_REGION_NAME = 'us-east-1' # MinIO 不需要区域，但 boto3 可能需要
#    AWS_S3_SIGNATURE_VERSION = 's3v4'
#    STATIC_URL = f'{AWS_S3_ENDPOINT_URL}/static/' # 如果使用 MinIO 存储静态文件
#    MEDIA_URL = f'{AWS_S3_ENDPOINT_URL}/media/'
#
#    **重要**: 确保 Django 容器可以访问到 PostgreSQL 和 MinIO。
#    如果 PostgreSQL 和 MinIO 也在 Docker 中运行，并且在同一个 Docker 网络中，
#    可以使用容器名作为 HOST (例如 'minio' 而不是 'YOUR_SERVER_IP')。
#    如果 PostgreSQL 在宿主机上运行，Django 容器可以使用宿主机的 IP 或 `host.docker.internal` (某些 Docker 版本)。
```

### **5.2 创建 Dockerfile**

在您的 Django 项目根目录下创建一个名为 `Dockerfile` 的文件：

```bash
# Dockerfile
FROM python:3.10-slim

# 设置工作目录
WORKDIR /app

# 设置 pip 国内镜像源 (可选, 加快构建速度)
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple

# 复制依赖文件并安装
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 复制项目代码到工作目录
COPY . .

# 暴露 Django 应用运行的端口
EXPOSE 8000

# 运行 Django 开发服务器 (生产环境推荐使用 Gunicorn 或 uWSGI)
# CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
# 使用 Gunicorn (确保 gunicorn 在 requirements.txt 中)
CMD ["gunicorn", "your_project_name.wsgi:application", "--bind", "0.0.0.0:8000"]
# 将 your_project_name 替换为您的 Django 项目的实际名称 (wsgi.py 所在的目录名)
```

### **5.3 构建并运行 Django Docker 镜像**

在包含 `Dockerfile` 的 Django 项目根目录下执行：

```bash
# 构建 Docker 镜像 (将 django-backend 替换为您的镜像名)
docker build -t django-backend .

# 运行 Django 容器
# 确保旧的同名容器已停止并移除:
# docker stop django-app && docker rm django-app
docker run -d \
  -p 8000:8000 \
  -e DJANGO_SETTINGS_MODULE=your_project_name.settings \
  --name django-app \
  django-backend
# 将 your_project_name.settings 替换为您的 Django settings 文件路径
```

## **第 6 步：配置 Nginx 反向代理**

Nginx 将作为前端服务器，接收外部请求并将其转发到 Django 应用。

```bash
# 1. 创建 Nginx 配置文件
# 将 YOUR_SERVER_IP 替换为您的服务器公网 IP 或域名
sudo bash -c "cat <<EOF > /etc/nginx/sites-available/django_proxy
server {
    listen 80;
    server_name YOUR_SERVER_IP; # 例如 123.45.6.78 或 example.com

    location / {
        proxy_pass http://127.0.0.1:8000; # Django 应用运行的地址和端口
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }

    location /static/ { # 如果 Django 自己处理静态文件 (DEBUG=True)
        alias /srv/django-app/staticfiles/; # 替换为你的 Django 项目静态文件收集目录
    }

    location /media/ { # 如果 Django 自己处理媒体文件 (DEBUG=True)
        alias /srv/django-app/media/; # 替换为你的 Django 项目媒体文件目录
    }
}
EOF"

# 2. 创建符号链接以启用该配置
# 先删除可能存在的默认配置的符号链接（如果它占用了 default_server）
# sudo rm /etc/nginx/sites-enabled/default
sudo ln -s /etc/nginx/sites-available/django_proxy /etc/nginx/sites-enabled/django_proxy

# 3. 测试 Nginx 配置
sudo nginx -t

# 4. 重启 Nginx 服务
sudo systemctl restart nginx
```

> **可视化操作与安全组说明 (Nginx & Django)**
>
> 1. **云服务器安全组**:
>    - 确保 HTTP 端口 `80/TCP` 已对公网开放。
>    - 如果未来配置 HTTPS，也需要开放 `443/TCP`。
> 2. **访问您的应用**:
>    - 在浏览器中输入 `http://YOUR_SERVER_IP` (例如 `http://123.45.6.78`)。
>    - 如果一切配置正确，您应该能看到您的 Django 应用首页。
>    - 如果部署了 Swagger，可以尝试访问 `http://YOUR_SERVER_IP:8000/api/swagger/` (路径取决于您的 Django URL 配置)。

## **第 7 步：自动化部署脚本**

将以下脚本保存为 `deploy_django_stack.sh`，赋予执行权限 (`chmod +x deploy_django_stack.sh`)，然后运行它。 **请务必在使用前仔细阅读脚本内容，并根据您的实际情况修改占位符和配置！**

```bash
#!/bin/bash

# --- 配置变量 (请务必根据您的实际情况修改!) ---
SERVER_IP="123.45.6.78" # 您的服务器公网 IP
POSTGRES_PASSWORD="YourSecurePostgresPassword!"
MINIO_ROOT_USER="admin"
MINIO_ROOT_PASSWORD="YourSecureMinioPassword!"
DJANGO_PROJECT_NAME="backend" # 您 Django 项目中包含 wsgi.py 的目录名
# DJANGO_REPO_URL="https://github.com/your-username/your-backend-repo.git" # 您的 Django 代码仓库地址
# DJANGO_PROJECT_DIR="/srv/django-app" # Django 项目将被克隆/放置到的目录

echo "--- 服务器部署脚本 ---"
echo "服务器 IP 将被设置为: $SERVER_IP"
echo "PostgreSQL 'postgres' 用户密码将设置为: $POSTGRES_PASSWORD"
echo "MinIO 管理员用户: $MINIO_ROOT_USER"
echo "MinIO 管理员密码: $MINIO_ROOT_PASSWORD"
echo "Django 项目名 (用于 Gunicorn): $DJANGO_PROJECT_NAME"
# echo "Django 代码仓库: $DJANGO_REPO_URL"
# echo "Django 项目目录: $DJANGO_PROJECT_DIR"
read -p "确认以上信息正确并开始部署吗？(yes/no): " confirmation
if [ "$confirmation" != "yes" ]; then
    echo "部署已取消。"
    exit 1
fi

echo "--- 1. 系统初始化与依赖安装 ---"
sudo apt update && sudo apt upgrade -y
sudo apt install -y ca-certificates curl gnupg lsb-release nginx git

echo "--- 2. 安装 Docker CE 和 Docker Compose ---"
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt update
sudo apt install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo mkdir -p /etc/docker
sudo tee /etc/docker/daemon.json <<EOF
{
  "registry-mirrors": [
    "https://docker.mirrors.ustc.edu.cn",
    "https://hub-mirror.c.163.com",
    "https://registry.docker-cn.com"
  ]
}
EOF
sudo systemctl daemon-reload
sudo systemctl restart docker
sudo usermod -aG docker $USER
echo "Docker 安装完成。请重新登录或执行 'newgrp docker' 以应用 docker 组权限。"
echo "按 Enter 继续..."
read

echo "--- 3. 安装 PostgreSQL 并配置 ---"
sudo apt install -y postgresql postgresql-contrib
PG_VERSION=$(pg_lsclusters | awk 'NR==2 {print $1}')
if [ -z "$PG_VERSION" ]; then
    echo "错误：无法检测到 PostgreSQL 版本。请手动配置。"
    exit 1
fi
echo "检测到 PostgreSQL 版本: $PG_VERSION"
sudo sed -i "s/#listen_addresses = 'localhost'/listen_addresses = '*'/g" /etc/postgresql/$PG_VERSION/main/postgresql.conf
sudo sh -c "echo 'host    all             all             0.0.0.0/0               md5' >> /etc/postgresql/$PG_VERSION/main/pg_hba.conf"
sudo systemctl restart postgresql
sudo -u postgres psql -c "ALTER USER postgres WITH PASSWORD '$POSTGRES_PASSWORD';"
echo "PostgreSQL 安装和配置完成。"

echo "--- 4. 部署 MinIO 对象存储 ---"
sudo mkdir -p /minio/data
sudo chmod -R 777 /minio/data # 确保 Docker 有权限写入
docker stop minio || true && docker rm minio || true # 确保旧容器被移除
docker run -d \
  --name minio \
  -p 9000:9000 \
  -p 9001:9001 \
  -v /minio/data:/data \
  -e "MINIO_ROOT_USER=${MINIO_ROOT_USER}" \
  -e "MINIO_ROOT_PASSWORD=${MINIO_ROOT_PASSWORD}" \
  quay.io/minio/minio:latest \
  server /data --console-address ":9001"
echo "MinIO 部署完成。请访问 http://$SERVER_IP:9001 并使用以下凭据登录："
echo "用户名: $MINIO_ROOT_USER"
echo "密码: $MINIO_ROOT_PASSWORD"
echo "登录后，请手动创建 'media' 和 'static' 存储桶，并根据需要设置其访问策略为公开可读。"
echo "按 Enter 继续..."
read

echo "--- 5. 部署 Django 后端应用 ---"
# echo "克隆 Django 项目代码从 $DJANGO_REPO_URL 到 $DJANGO_PROJECT_DIR..."
# sudo mkdir -p $DJANGO_PROJECT_DIR
# sudo chown $USER:$USER $DJANGO_PROJECT_DIR # 确保当前用户有权限
# git clone $DJANGO_REPO_URL $DJANGO_PROJECT_DIR
# cd $DJANGO_PROJECT_DIR
echo "假设 Django 项目代码已位于当前目录或指定目录。"
echo "请确保您的 Django 项目 (例如: ./backend/settings.py) 已配置好数据库和 MinIO 连接。"
echo "数据库主机应指向 '$SERVER_IP' (如果 PostgreSQL 在宿主机上运行)。"
echo "MinIO Endpoint URL 应指向 'http://$SERVER_IP:9000'。"
echo "按 Enter 继续以创建 Dockerfile 并构建镜像..."
read

# 创建 Dockerfile
cat <<EOF > Dockerfile
FROM python:3.10-slim
ENV PYTHONUNBUFFERED 1
WORKDIR /app
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
EXPOSE 8000
# 将 your_project_name 替换为您的 Django 项目名 (wsgi.py 所在的目录名)
CMD ["gunicorn", "${DJANGO_PROJECT_NAME}.wsgi:application", "--bind", "0.0.0.0:8000"]
EOF
echo "Dockerfile 已创建。"

echo "构建 Django Docker 镜像 (django-backend)..."
docker build -t django-backend .

echo "运行 Django Docker 容器 (django-app)..."
docker stop django-app || true && docker rm django-app || true # 确保旧容器被移除
docker run -d \
  -p 8000:8000 \
  -e DJANGO_SETTINGS_MODULE=${DJANGO_PROJECT_NAME}.settings \
  --name django-app \
  django-backend
echo "Django 应用容器已启动。"
echo "按 Enter 继续配置 Nginx..."
read

echo "--- 6. 配置 Nginx 反向代理 ---"
# Nginx 已在步骤1安装
sudo bash -c "cat <<EOF > /etc/nginx/sites-available/django_proxy
server {
    listen 80;
    server_name $SERVER_IP;

    client_max_body_size 100M; # 允许上传大文件

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;
    }

    # 如果您希望 Nginx 直接处理静态文件和媒体文件 (生产环境推荐)
    # 请确保 Django 的 collectstatic 已将文件收集到 Nginx 可访问的路径
    # location /static/ {
    #     alias /path/to/your/django_project/staticfiles/; # 替换为实际路径
    # }
    # location /media/ {
    #     alias /path/to/your/django_project/mediafiles/; # 替换为实际路径
    # }
}
EOF"
# 确保旧的 default 站点（如果存在且冲突）被禁用
if [ -L /etc/nginx/sites-enabled/default ]; then
    sudo rm /etc/nginx/sites-enabled/default
fi
# 强制创建或更新符号链接
sudo ln -sf /etc/nginx/sites-available/django_proxy /etc/nginx/sites-enabled/django_proxy

echo "测试 Nginx 配置..."
sudo nginx -t
if [ \$? -ne 0 ]; then
    echo "Nginx 配置测试失败！请检查 /etc/nginx/sites-available/django_proxy 文件。"
    exit 1
fi
echo "重启 Nginx 服务..."
sudo systemctl restart nginx
echo "Nginx 配置完成。"

echo "--- 部署完成！---"
echo "请确保您的云服务器安全组已开放以下端口："
echo "  - PostgreSQL: 5432/TCP (如果需要远程访问数据库)"
echo "  - MinIO API: 9000/TCP"
echo "  - MinIO 控制台: 9001/TCP"
echo "  - HTTP (Nginx): 80/TCP"
echo ""
echo "您现在应该可以通过 http://$SERVER_IP 访问您的 Django 应用。"
echo "MinIO 控制台: http://$SERVER_IP:9001"
echo "如果 Django 应用需要数据库迁移，请在容器内执行:"
echo "  docker exec -it django-app python manage.py makemigrations your_app_name"
echo "  docker exec -it django-app python manage.py migrate"
echo "如果需要创建超级用户:"
echo "  docker exec -it django-app python manage.py createsuperuser"
```

## **8. 用户数据迁移策略 (可选)**

当系统用户量增长，或者需要从旧系统迁移数据时，需要考虑数据迁移策略。

1. **新用户注册与数据处理**:
   - **Django Auth**: Django 自带的 `django.contrib.auth` 系统能很好地处理新用户注册、密码哈希存储、登录认证等。当新用户通过您的 API (例如 `/api/register/`) 注册时，会在 `auth_user` 表（或您自定义的用户模型表）中创建一条记录。
   - **关联用户信息 (`UserInfo`)**: 如果您有一个单独的 `UserInfo` 模型通过 `OneToOneField` 或 `ForeignKey` 关联到主用户模型，确保在用户注册成功后，或用户首次编辑个人资料时，创建或更新对应的 `UserInfo` 记录。`user_id` 将作为关联两个表的键。
2. **从旧系统迁移数据 (如果适用)**:
   - **数据导出**: 从旧系统导出用户数据，通常为 CSV、JSON 或 SQL dump 格式。
   - **数据清洗与转换**: 清洗数据，使其符合新系统的数据模型。特别注意密码的处理，如果旧系统密码哈希算法与 Django 不兼容，用户可能需要在首次登录新系统时重置密码。
   - **数据导入**:
     - **Django 管理命令**: 编写自定义的 Django management command (`python manage.py your_custom_command`)，使用 Django ORM 来创建用户和关联信息。这是推荐的方式，因为它会处理所有模型逻辑和信号。
     - **直接 SQL**: 对于非常大的数据集，直接使用 SQL 导入到 PostgreSQL 可能更快，但需要非常小心，确保数据完整性和关联正确，并且后续可能需要手动处理 Django 的 contenttypes 等。
     - **第三方库**: 例如 `django-import-export` 库可以帮助处理复杂的数据导入导出。
3. **数据完整性与验证**:
   - 在导入过程中，使用 Django 模型的 `full_clean()` 方法或表单验证来确保数据符合新模型的约束。
   - 利用数据库的约束（如 `UNIQUE`, `NOT NULL`）来保证数据质量。
4. **处理静态文件和媒体文件 (如头像)**:
   - 如果旧系统有用户上传的文件，需要将这些文件迁移到新的存储位置（例如 MinIO）。
   - 更新数据库中指向这些文件的路径或 URL。
   - 如果文件名或路径结构发生变化，需要编写脚本来批量更新。
5. **扩展性考虑**:
   - **数据库**: PostgreSQL 本身具有良好的扩展性。对于非常大的负载，可以考虑读写分离、分区等策略。
   - **Django 应用**: 使用 Gunicorn 或 uWSGI 配合多个 worker 进程可以处理更多并发请求。Nginx 作为反向代理可以进行负载均衡。
   - **缓存**: 对常用数据和计算结果使用缓存（如 Redis、Memcached）可以显著提高性能。
6. **备份与恢复**:
   - **数据库**: 定期使用 `pg_dump` 备份 PostgreSQL 数据库。制定恢复计划。
   - **MinIO**: 定期备份 MinIO 的数据卷 (`/minio/data`)。MinIO 也支持自身的复制和纠删码功能来提高数据可靠性。
   - **应用代码和配置**: 使用版本控制 (如 Git) 管理代码，并备份 Docker 镜像和相关配置文件。

## **9. MinIO 和 Django (Docker) Debug 指南**

在部署和运行 MinIO 及 Docker化的 Django 应用时，可能会遇到一些常见问题。以下是一些调试步骤和技巧，结合了您之前遇到的情况：

### **9.1 MinIO Client (`mc`) 相关问题**

1. **`mc`: Segmentation fault**:

   - **原因**: `mc` 二进制文件损坏、与系统架构不兼容（例如在 ARM 服务器上运行了 AMD64 版本）或下载不完整。

   - **解决方案**:

     1. **彻底卸载旧/损坏的 `mc`**:

        ```bash
        sudo rm -f /usr/local/bin/mc
        which mc # 确认已删除，应无输出
        ```

     2. **下载正确的官方版本** (假设为 linux-amd64):

        ```bash
        wget https://dl.min.io/client/mc/release/linux-amd64/mc
        chmod +x mc
        sudo mv mc /usr/local/bin/
        ```

     3. **验证**: `mc --version`

   - **检查文件类型**: `file /usr/local/bin/mc` (应显示 `ELF 64-bit LSB executable, x86-64,...`)

2. **`mc alias set ...`: The request signature we calculated does not match...**:

   - **原因**:
     - **Access Key / Secret Key 错误**: 您提供的 `admin` 和 `SecurePassword123!` (或您实际使用的密码) 与 MinIO 服务启动时配置的 `MINIO_ROOT_USER` / `MINIO_ROOT_PASSWORD` 不匹配。
     - **Endpoint URL 错误**: URL 格式不正确 (例如多了斜杠 `http://123.45.6.78/:9000`) 或地址/端口错误。正确应为 `http://YOUR_SERVER_IP:9000`。
     - MinIO 服务未运行或端口 9000 未正确映射/防火墙未开放。
   - **解决方案**:
     1. 确认 MinIO 容器正在运行且端口 9000 已映射: `docker ps | grep minio`
     2. 仔细核对 `mc alias set` 命令中的 Access Key (用户名), Secret Key (密码), 和 Endpoint URL。
     3. 确保密码中没有特殊字符导致命令行解析问题，或者用单引号包裹密码：`mc alias set myminio http://123.45.6.78:9000 admin 'YourSecureMinioPassword!'`

3. **`mc admin user password ...`: `password` is not a recognized command**:

   - **原因**:
     - `mc` 版本过旧，不支持该子命令。
     - 命令语法错误。
   - **解决方案**:
     1. **升级 `mc`**: 确保您使用的是最新版本的 `mc` (参考上面安装步骤)。
     2. **正确语法**: `mc admin user password ALIAS USERNAME NEW_PASSWORD` 例如: `mc admin user password myminio admin 'NewSecurePassword123!'` (确保 `myminio` 别名已成功设置)。

### **9.2 Docker 相关问题**

1. **`docker run ...`: address already in use (e.g., for port 8000)**:

   - **原因**: 您尝试映射到宿主机的端口 (例如 8000) 已经被其他进程占用。

   - **解决方案**:

     1. 查找并停止占用端口的进程:

        ```bash
        sudo lsof -i :8000
        # 或
        sudo netstat -tulnp | grep 8000
        # 找到 PID 后，使用 sudo kill <PID> 或 sudo kill -9 <PID>
        ```

     2. 如果占用者是另一个 Docker 容器，先停止并移除它:

        ```bash
        docker ps -a # 查看所有容器，找到占用端口的容器
        docker stop <container_id_or_name>
        docker rm <container_id_or_name>
        ```

     3. 或者，更改您当前要运行的容器的端口映射，例如将 Django 映射到宿主机的 8001 端口： `docker run -p 8001:8000 ...` (同时需要更新 Nginx 配置中的 `proxy_pass` 指向 `http://127.0.0.1:8001;`)

2. **`docker run ...`: Conflict. The container name "/django-app" is already in use...**:

   - **原因**: 已存在一个同名的 Docker 容器 (即使它已停止)。
   - **解决方案**:
     1. 删除旧的同名容器: `docker rm django-app` (如果容器已停止) 或 `docker stop django-app && docker rm django-app` (如果正在运行)。
     2. 或者，为新容器指定一个不同的名称：`docker run --name django-app-v2 ...`

3. **`docker run ...`: invalid reference format`或`django-backend: command not found`**:

   - **原因**: Docker 无法找到您指定的镜像 `django-backend`。
     - 镜像名称拼写错误或大小写不匹配 (Docker 镜像名通常是小写)。
     - 镜像尚未成功构建，或者构建时使用了不同的标签。
   - **解决方案**:
     1. 确认镜像是否存在且名称正确: `docker images | grep django-backend`
     2. 如果不存在，请确保在 Django 项目根目录 (包含 `Dockerfile`) 下重新构建: `docker build -t django-backend .`
     3. 确保 `docker run` 命令中使用的镜像名称与 `docker images` 中显示的完全一致。

### **9.3 Django (Docker 内部) 调试**

1. **检查 Django 容器状态和日志**:

   ```bash
   docker ps -a | grep django-app  # 查看容器是否正在运行
   docker logs django-app          # 查看 Django 容器的实时日志（非常重要！）
   ```

   - 日志会显示 Gunicorn/Django 的启动信息、任何 Python 错误、数据库连接问题等。

2. **进入 Django 容器内部进行调试**:

   ```bash
   docker exec -it django-app bash # 进入容器的 shell
   ```

   - 进入容器后，您可以：
     - 检查文件是否存在，路径是否正确。
     - 手动运行 `python manage.py check` 查看是否有配置问题。
     - 尝试连接数据库：`python manage.py dbshell` (如果安装了 psql 客户端)。
     - 查看环境变量是否正确设置。

3. **数据库连接问题**:

   - **错误**: 日志中可能出现 `OperationalError: could not connect to server` 或类似错误。
   - **检查**:
     - PostgreSQL 服务是否在宿主机或另一容器中运行。
     - Django `settings.py` 中的 `DATABASES` 配置（`HOST`, `PORT`, `NAME`, `USER`, `PASSWORD`）是否正确。
     - 如果 PostgreSQL 在宿主机，Django 容器是否能访问到宿主机的 IP 和端口。可能需要在 `docker run` 时使用 `--network="host"` (不推荐，除非必要) 或确保 Docker 网络配置允许。更常见的是使用服务器的公网 IP (确保 PostgreSQL 监听 `*` 且防火墙允许)。
     - PostgreSQL 的 `pg_hba.conf` 是否允许来自 Docker 容器 IP 地址范围的连接。

### **9.4 Nginx 调试**

1. **`nginx: [emerg] invalid number of arguments in "proxy_set_header" directive...`**:
   - **原因**: `proxy_set_header` 指令语法错误，通常是参数数量不对或变量名错误。
   - **示例错误**: `proxy_set_header X-Forwarded-For $remote_addr;`
   - **正确示例**: `proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;`
   - **解决方案**: 仔细检查 Nginx 配置文件 (`/etc/nginx/sites-available/django_proxy`) 中的所有 `proxy_set_header` 指令，确保它们都遵循 `proxy_set_header <字段名> <值>;` 的格式。
   - **注意变量**: Nginx 内置变量如 `$host`, `$remote_addr`, `$proxy_add_x_forwarded_for` 不需要额外的转义符。
2. **Nginx 502 Bad Gateway**:
   - **原因**: Nginx 无法连接到 `proxy_pass` 指令中指定的后端 Django 应用。
   - **解决方案**:
     1. **确认 Django 容器运行**: `docker ps | grep django-app`。
     2. **检查 Django 容器日志**: `docker logs django-app`，查看是否有启动错误。
     3. **确认 Django 监听端口**: 确保 Django 应用 (Gunicorn) 在容器内部监听 `0.0.0.0:8000`。
     4. **确认 Nginx `proxy_pass` 地址**: 通常是 `http://127.0.0.1:8000;` (因为 Django 容器的 8000 端口映射到了宿主机的 8000 端口，Nginx 从宿主机访问此端口)。
     5. **网络测试**: 在服务器上尝试 `curl http://127.0.0.1:8000`，看是否能访问到 Django 应用。
3. **Django 404 for `/swagger/` (但 `/api/swagger/` works)**:
   - **原因**: 访问的 URL 路径与 Django `urls.py` 中定义的路由不匹配。
   - **解决方案**: 确保您在浏览器中访问的是 Django `urls.py` 中为 Swagger UI 定义的正确路径，例如 `http://YOUR_SERVER_IP/api/swagger/`。

### **9.5 通用调试技巧**

- **逐步验证**: 从底层服务（PostgreSQL, MinIO）开始，确保它们独立运行时正常，然后再验证 Django 应用，最后是 Nginx。

- **简化配置**: 如果遇到复杂问题，尝试暂时简化配置（例如，移除 Nginx，直接暴露 Django 容器端口进行测试）以缩小问题范围。

- **查看所有日志**: 同时关注 PostgreSQL, MinIO, Django (Gunicorn), Nginx 的日志，它们通常会提供关键的错误信息。

- **网络工具**: 使用 `ping`, `curl`, `telnet`, `netstat`, `ss` 等工具检查网络连通性和端口监听情况。

  ```bash
  # 检查端口是否被监听
  sudo netstat -tulnp | grep 5432 # PostgreSQL
  sudo netstat -tulnp | grep 9000 # MinIO API
  sudo netstat -tulnp | grep 9001 # MinIO Console
  sudo netstat -tulnp | grep 8000 # Django App / Nginx
  sudo netstat -tulnp | grep 80   # Nginx
  ```

通过这些步骤，您应该能够定位并解决部署过程中遇到的大部分问题。