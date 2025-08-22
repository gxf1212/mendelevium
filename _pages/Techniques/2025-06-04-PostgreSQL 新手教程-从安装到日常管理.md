---
title: "PostgreSQL Beginner Tutorial: From Installation to Daily Management"
date: "2025-06-04"
tags: [postgresql, database, tutorial, installation, administration, sql, ubuntu]
---


## **PostgreSQL 新手教程：从安装到日常管理**

本教程旨在为初学者提供一个清晰、易懂的 PostgreSQL 上手指南，内容涵盖安装、初始化、用户创建、数据库操作以及常见问题解决。

### **1. 📦 安装 PostgreSQL**

在开始之前，建议先更新您的系统包列表。以下命令适用于基于 RHEL/CentOS 的系统（使用 `yum`）：

```
sudo yum update -y
sudo yum install postgresql postgresql-server postgresql-contrib -y
```

- `postgresql`: 包含客户端程序。
- `postgresql-server`: 包含 PostgreSQL 服务器本身。
- `postgresql-contrib`: 包含一些额外的贡献模块和工具。

安装过程中，如果系统提示有关依赖项的更新，请确认并继续安装。完成后，您可能会看到类似以下的输出，表明相关软件包已成功安装或更新：

```
Dependency Updated:
  postgresql-devel.x86_64 0:9.2.24-9.el7_9
  postgresql-docs.x86_64 0:9.2.24-9.el7_9
  postgresql-libs.x86_64 0:9.2.24-9.el7_9
  postgresql-server.x86_64 0:9.2.24-9.el7_9
Complete!
```

**注意**：对于其他 Linux 发行版（如 Debian/Ubuntu），安装命令会有所不同（例如 `sudo apt-get install postgresql postgresql-contrib`）。

### **2. 🚀 初始化数据库集群**

这是**非常关键的一步**！PostgreSQL 安装完成后，并不会自动创建一个可用的数据库环境。您需要首先初始化一个“数据库集群”。一个数据库集群是一组数据库的集合，由单个 PostgreSQL 服务器实例管理。

```
sudo postgresql-setup initdb
```

- 此命令会在默认位置（通常是 `/var/lib/pgsql/data`）创建数据库集群所需的数据目录结构。
- 它还会生成一些核心的配置文件，例如 `postgresql.conf`（主配置文件，控制服务器行为）和 `pg_hba.conf`（客户端认证配置文件，控制哪些用户可以从哪些主机连接以及如何认证）。
- 如果这个步骤被跳过或者数据目录未正确初始化，后续启动 PostgreSQL 服务将会失败。

### **3. ▶️ 启动 PostgreSQL 服务**

初始化数据库集群后，就可以启动 PostgreSQL 服务了：

```
sudo systemctl start postgresql
```

为了让 PostgreSQL 在系统重启后自动启动（推荐做法），可以设置开机自启：

```
sudo systemctl enable postgresql
```

### **4. 🩺 验证服务状态**

启动服务后，检查其运行状态以确保一切正常：

```
systemctl status postgresql.service
```

如果服务正常运行，您应该在输出中看到 `Active: active (running)` 字样以及服务的启动时间等信息。

### **5. 🔑 连接到 PostgreSQL (理解登录机制)**

当您首次安装并初始化 PostgreSQL 后，系统会自动创建一个名为 `postgres` 的**操作系统用户**和一个同名的 **PostgreSQL 超级用户（数据库角色）**。

**默认的登录方式（Peer Authentication）：**

在很多默认配置中，PostgreSQL 使用“peer”认证方式进行本地连接。这意味着，如果您的操作系统用户名与您尝试连接的 PostgreSQL 用户名相同，并且您有权访问数据库套接字文件，那么 PostgreSQL 会信任操作系统已经验证了您的身份，从而允许您无需密码直接登录。

这就是为什么我们使用以下命令以 `postgres` 操作系统用户的身份来连接 `postgres` 数据库超级用户：

```
sudo -u postgres psql
```

- `sudo -u postgres`: 这部分命令的含义是“以操作系统用户 `postgres` 的身份执行后续命令”。
- `psql`: 这是 PostgreSQL 的命令行交互式客户端工具。

执行此命令后，您可能会看到类似以下的提示：

```
could not change directory to "/root"
psql (9.2.24)
Type "help" for help.

postgres=#
```

- `could not change directory to "/root"`: 这个提示通常可以忽略。它是因为您通过 `sudo` 切换到 `postgres` 用户，但 `postgres` 用户可能没有权限访问您当前所在的目录（例如 `/root`）。这不影响 `psql` 的连接和操作。
- `postgres=#`: 这个提示符表明您已成功以 `postgres` 超级用户身份连接到了 PostgreSQL。`postgres=` 左边的 `postgres` 通常表示您当前连接的默认数据库名（也叫 `postgres`），`=` 表示您是超级用户 (`#` 也是超级用户的标志，不同版本或配置可能略有不同)。

**关于指定数据库登录：**

当您运行 `psql` 时，如果您不指定要连接的数据库，它会尝试连接一个与您当前操作系统用户名同名的数据库。由于我们使用了 `sudo -u postgres psql`，它会尝试以数据库用户 `postgres` 连接名为 `postgres` 的数据库。这个 `postgres` 数据库是一个在初始化时自动创建的默认管理用数据库。

**简单来说，用户登录数据库时，PostgreSQL 会根据 `pg_hba.conf` 文件中的规则来验证您的身份。** `psql` 命令本身可以带有参数来指定用户名、数据库名、主机等，例如：

```
psql -U <数据库用户名> -d <数据库名> -h <主机地址>
```

如果您不写这些参数，`psql` 会使用一些默认值（通常是当前操作系统用户名作为数据库用户名和数据库名）。

### **6. ✨ `postgres` 超级用户的特权和功能**

在 PostgreSQL 中，`postgres` 用户（或任何被赋予 `SUPERUSER` 属性的角色）拥有数据库集群内的最高权限。这些特权和功能包括但不限于：

- **不受限制的访问权限**：可以访问集群中的所有数据库和所有对象（表、视图、函数等），无视常规的权限检查。
- **创建和管理数据库**：可以创建新的数据库 (`CREATEDB` 属性的体现)。
- **创建和管理角色（用户和组）**：可以创建、修改、删除其他用户和用户组 (`CREATEROLE` 属性的体现)。
- **修改配置参数**：可以更改 PostgreSQL 服务器的配置参数（`postgresql.conf` 中的设置）。
- **执行特权操作**：例如加载C语言函数、执行某些诊断或维护命令、绕过某些安全限制等。
- **复制权限**：可以启动和管理流复制 (`REPLICATION` 属性的体现)。
- **绕过行级安全策略 (Row-Level Security)**。

**正是因为超级用户权限过大，通常不建议在日常应用程序中使用超级用户账户。** 最佳实践是为应用程序创建具有完成其任务所需最小权限的普通用户。

### **7. 👤 创建新用户**

在 `psql` 提示符 (`postgres=#`) 下，您可以使用 `postgres` 超级用户的权限来创建新的数据库用户（在 PostgreSQL 中称为“角色”）。

```
CREATE USER wangpeng WITH CREATEDB PASSWORD 'your_strong_password';
```

- `CREATE USER wangpeng`: 创建一个名为 `wangpeng` 的新用户。
- `WITH CREATEDB`: 授予该用户创建新数据库的权限。这是一个常见的权限，但不等同于超级用户。其他可选权限包括 `SUPERUSER`, `CREATEROLE`, `LOGIN` (默认就有), `REPLICATION` 等。
- `PASSWORD 'your_strong_password'`: 为新用户设置一个密码。**请务必将其替换为一个强密码。** 如果不指定 `PASSWORD`，用户将无法通过密码认证登录（可能需要其他认证方式）。

### **8. 📋 查看用户列表**

要查看当前数据库集群中存在的所有角色（用户和组）及其属性，可以在 `psql` 中使用 `\du` 元命令：

```
\du
```

输出示例：

```
                                     List of roles
 Role name |                   Attributes                   | Member of 
-----------+------------------------------------------------+-----------
 postgres  | Superuser, Create role, Create DB, Replication | {}
 wangpeng  | Create DB                                      | {}
```

- `Role name`: 用户名。
- `Attributes`: 该用户拥有的权限属性。
- `Member of`: 该用户所属的用户组（这里为空）。

### **9. 🚪 退出 psql**

当您完成在 `psql` 中的操作后，可以使用 `\q` 元命令退出：

```
\q
```

或者直接按 `Ctrl+D`。

### **10. 🔑 普通用户修改自己的密码**

一个已经创建的、具有登录权限的普通用户（例如我们前面创建的 `wangpeng`），可以在连接到数据库后修改自己的密码。

首先，该用户需要能够登录。假设 `pg_hba.conf` 已配置为允许密码认证（例如 `md5`），用户 `wangpeng` 可以这样登录（可能需要先退出 `postgres` 用户的 `psql` 会话）：

```
psql -U wangpeng -d postgres -h localhost -W
```

(这里连接到 `postgres` 数据库只是为了执行 `ALTER USER` 命令，也可以连接到该用户有权限的其他数据库)

然后，在 `psql` 提示符下（此时应该是 `wangpeng=>` 或类似），执行：

```
ALTER USER wangpeng WITH PASSWORD 'new_strong_password';
```

这样，用户 `wangpeng` 就成功修改了自己的密码。

**注意**：普通用户只能修改自己的密码，不能修改其他用户的密码。只有超级用户或具有 `CREATEROLE` 权限的用户才能修改其他用户的密码或属性。

## **PostgreSQL 新手教程：创建和管理数据库**

### **1. 创建数据库**

一旦您以具有 `CREATEDB` 权限的用户（例如 `postgres` 超级用户或我们之前创建的 `wangpeng`）登录到 `psql`，就可以创建新的数据库了：

```
-- 假设以 wangpeng 用户登录后执行
CREATE DATABASE mydatabase;
```

### **2. 为数据库创建专属普通用户 (推荐)**

通常，为每个应用程序或主要功能创建一个专用的数据库用户是一个好习惯，而不是直接使用像 `wangpeng` 这样可能具有创建数据库权限的用户。

```
-- 假设仍以 postgres 或 wangpeng (有 CREATEROLE 潜质，或 postgres 执行) 登录
CREATE USER myapp_user WITH PASSWORD 'another_strong_password';
```

这里我们创建了一个名为 `myapp_user` 的用户，它默认只具有 `LOGIN` 权限。

### **3. 授予用户对特定数据库的权限**

新创建的 `myapp_user` 默认情况下对我们刚创建的 `mydatabase` 没有任何操作权限（除了连接，如果认证方式允许）。我们需要明确授予它权限：

```
-- 授予 myapp_user 对 mydatabase 数据库的所有基本权限
GRANT ALL PRIVILEGES ON DATABASE mydatabase TO myapp_user;

-- (可选) 如果希望 myapp_user 能够在该数据库中创建表等对象，
-- 可能还需要更改该数据库中默认 schema (如 public) 的权限，
-- 或者为该用户创建一个专属的 schema 并授予权限。
-- 例如，允许在新数据库的 public schema 中创建对象：
-- \c mydatabase  -- 首先连接到目标数据库
-- GRANT CREATE ON SCHEMA public TO myapp_user;
-- GRANT USAGE ON SCHEMA public TO myapp_user; -- 允许使用 public schema
```

- `GRANT ALL PRIVILEGES ON DATABASE mydatabase TO myapp_user;`: 这授予了用户在 `mydatabase` 上的连接权限，以及在该数据库内创建 schema 的权限（如果默认权限允许）。但**不包括**在该数据库内创建表、序列等对象的权限，也不包括对已有对象的访问权限。对具体对象（表、视图、序列、函数等）的权限需要单独授予 (例如 `GRANT SELECT, INSERT ON my_table TO myapp_user;`)。

### **4. 使用新用户登录特定数据库**

现在，`myapp_user` 可以尝试登录 `mydatabase` 了：

```
# 在操作系统的命令行中执行
psql -U myapp_user -d mydatabase -h localhost -W
```

- `-U myapp_user`: 指定要使用的数据库用户名。
- `-d mydatabase`: 指定要连接的数据库名。
- `-h localhost`: 指定数据库服务器的主机地址（如果是本地服务器，此参数有时可以省略，取决于配置）。
- `-W`: 强制提示输入密码。

成功登录后，提示符会变为类似 `mydatabase=>` （如果 `myapp_user` 不是超级用户）。

## **完整流程示例（概览）**

以下是一个简化的完整流程，展示了从首次安装后的一系列操作：

```
# (在操作系统命令行中)
# 1. 初始化数据库（仅限首次安装 PostgreSQL 后执行一次）
sudo postgresql-setup initdb

# 2. 启动 PostgreSQL 服务
sudo systemctl start postgresql
sudo systemctl enable postgresql # 设置开机自启

# 3. 切换到 postgres 操作系统用户并进入 psql
sudo -u postgres psql

# --- 以下命令在 psql (postgres=#) 提示符下执行 ---
# 4. 创建一个新角色（用户）并赋予创建数据库的权限
CREATE USER db_admin WITH CREATEDB PASSWORD 'admin_password123';

# 5. (可选) db_admin 用户退出并重新登录，或继续以 postgres 操作
-- \q  -- 退出 postgres 的 psql 会话

-- (如果退出了，则重新以 db_admin 登录)
-- psql -U db_admin -d postgres -W  (输入 admin_password123)

-- 6. 创建一个新的数据库 (假设由 db_admin 创建)
CREATE DATABASE new_application_db;

# 7. 为应用程序创建一个权限更受限的用户
CREATE USER app_user WITH PASSWORD 'app_password456';

# 8. 授予 app_user 连接到新数据库的权限
GRANT CONNECT ON DATABASE new_application_db TO app_user;

# 9. 切换到新数据库以授予更细致的权限
\c new_application_db

# 10. 授予 app_user 在 public schema 中创建对象和使用 schema 的权限
GRANT CREATE ON SCHEMA public TO app_user;
GRANT USAGE ON SCHEMA public TO app_user;
-- (如果需要，还可以授予对特定表的 SELECT, INSERT, UPDATE, DELETE 权限)
-- 例如: GRANT SELECT, INSERT ON TABLE my_table TO app_user;

# 11. 退出 psql
\q
# --- psql 操作结束 ---

# (可选) 12. 修改认证方式以允许密码登录 (如果默认是 peer)
# 编辑 /var/lib/pgsql/data/pg_hba.conf 文件
# sudo vi /var/lib/pgsql/data/pg_hba.conf
# 将相关 'local' 或 'host' 行的认证方法从 'peer' 或 'ident' 改为 'md5' 或 'scram-sha-256'
# 例如:
# local   all             all                                     md5
# host    all             all             127.0.0.1/32            md5
# host    all             all             ::1/128                 md5

# (如果修改了 pg_hba.conf) 13. 重启 PostgreSQL 服务使配置生效
# sudo systemctl restart postgresql

# 14. 以新创建的 app_user 登录到其专属数据库
# psql -U app_user -d new_application_db -h localhost -W
# (输入 app_password456)
```

## **常见问题及解决方案**

### **问题：PostgreSQL 服务启动失败**

- **现象**：执行 `sudo systemctl start postgresql` 时提示失败，日志中可能出现 `Failed at step EXEC_START pre spawning, Unit entered failed state.`
- 主要原因
  - **数据目录未初始化**：`/var/lib/pgsql/data` (或您系统上的默认数据目录) 不存在或为空。这是最常见的原因。
  - 服务配置错误或权限问题。
- 解决方法
  1. **确保已初始化数据库集群**：`sudo postgresql-setup initdb`
  2. 再次尝试启动服务：`sudo systemctl start postgresql`
  3. 检查服务状态和日志获取更详细错误信息：`systemctl status postgresql.service` 和 `journalctl -xeu postgresql.service`

### **问题：无法连接到 PostgreSQL (`psql: could not connect to server: No such file or directory`)**

- **现象**：执行 `psql` 时提示上述错误。
- 主要原因
  - PostgreSQL 服务未启动。
  - `psql` 尝试连接的 Unix 域套接字文件不存在（通常因为服务未运行或配置错误）。
- 解决方法
  1. **确保服务已启动**：`sudo systemctl status postgresql`，如果未运行则 `sudo systemctl start postgresql`。
  2. 检查数据目录是否已正确初始化（应包含 `PG_VERSION`, `global`, `pg_hba.conf` 等文件）：`ls /var/lib/pgsql/data`。

### **问题：认证失败（例如 `Peer authentication failed for user "username"`）**

- **现象**：尝试使用特定用户和密码登录时，即使密码正确，也提示认证失败。

- **主要原因**：`pg_hba.conf` 文件中配置的认证方法不允许您尝试的连接类型或用户使用密码认证。例如，对于本地连接，默认可能配置为 `peer` 或 `ident` 认证，而不是 `md5` 或 `scram-sha-256` (密码认证)。

- 解决方法

  修改 `pg_hba.conf` 文件

  ```
  sudo vi /var/lib/pgsql/data/pg_hba.conf 
  ```

  找到与您的连接类型（

  ```
  local
  ```

   表示 Unix 域套接字连接，

  ```
  host
  ```

   表示 TCP/IP 连接）、数据库、用户匹配的行，将其最后的认证方法修改为 

  ```
  md5
  ```

   (较旧，但兼容性好) 或 

  ```
  scram-sha-256
  ```

   (更安全，推荐用于新版本)。 例如，要允许所有本地用户通过密码连接所有数据库：

  ```
  # TYPE  DATABASE        USER            ADDRESS                 METHOD
  local   all             all                                     md5 
  host    all             all             127.0.0.1/32            md5
  host    all             all             ::1/128                 md5
  ```

  注意：修改 

  ```
  pg_hba.conf
  ```

   时要小心，不正确的配置可能导致无法连接。

  ```
  all all
  ```

  是一种比较宽松的配置，生产环境应根据需要进行限制。

  重启 PostgreSQL 服务使配置生效：

  ```
  sudo systemctl restart postgresql
  ```

### **问题：登录后提示 “FATAL: database "username" does not exist”**

- **现象**：使用某个用户（例如 `newuser`）通过 `psql -U newuser` 登录时，如果该用户对应的同名数据库 (`newuser`) 不存在，会报此错误。

- **原因**：`psql` 在未指定 `-d <数据库名>` 时，默认尝试连接与用户名同名的数据库。

- 解决方法

  1. 在登录时明确指定一个存在的数据库：`psql -U newuser -d postgres` (连接到默认的 `postgres` 数据库) 或 `psql -U newuser -d mydatabase` (连接到您已创建的 `mydatabase`)。

  2. 或者，如果您确实需要一个与该用户同名的数据库，请先以有权限的用户，如 

     ```
     postgres
     ```

     登录并创建它：

     ```
     -- 以 postgres 用户或其他有权限用户执行
     CREATE DATABASE newuser OWNER newuser; -- 创建数据库并将所有权赋给 newuser
     ```

## **附录：常用 PostgreSQL (psql) 命令速查**

| **操作**                 | **PostgreSQL (psql) 命令/ 命令**                             |
| ------------------------ | ------------------------------------------------------------ |
| 安装 PostgreSQL          | `sudo yum install postgresql postgresql-server postgresql-contrib -y` |
| 初始化数据库集群         | `sudo postgresql-setup initdb`                               |
| 启动服务                 | `sudo systemctl start postgresql`                            |
| 查看服务状态             | `systemctl status postgresql.service`                        |
| 以`postgres`用户进入psql | `sudo -u postgres psql`                                      |
| 列出所有数据库           | `\l` 或 `\list`                                              |
| 连接到特定数据库         | `\c database_name` 或 `\connect database_name`               |
| 列出当前数据库中的表     | `\dt`                                                        |
| 查看表结构               | `\d table_name`                                              |
| 列出所有用户(角色)       | `\du`                                                        |
| 列出所有schema           | `\dn`                                                        |
| 显示当前连接信息         | `\conninfo`                                                  |
| 执行上一条命令           | `\g`                                                         |
| 编辑上一条命令           | `\e`                                                         |
| 显示帮助信息             | `\?` (psql元命令帮助) 或 `help;` (SQL命令帮助)               |
| 退出 psql                | `\q`                                                         |
| **SQL 命令示例**         |                                                              |
| 创建用户                 | `CREATE USER username WITH PASSWORD 'password';`             |
| 修改用户密码             | `ALTER USER username WITH PASSWORD 'new_password';`          |
| 创建数据库               | `CREATE DATABASE dbname;`                                    |
| 删除数据库               | `DROP DATABASE dbname;`                                      |
| 授予权限                 | `GRANT ALL PRIVILEGES ON DATABASE dbname TO username;`       |
| 授予表权限               | `GRANT SELECT, INSERT, UPDATE, DELETE ON TABLE tablename TO username;` |
