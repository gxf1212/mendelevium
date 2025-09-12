# Commands

## jekyll
```shell
bundle install
bundle exec jekyll build
```

## posts

```shell
sed -i 's/\\\\/\\/g' *.md
sed -i 's/\\\_/\_/g' *.md
sed -i 's/\\\*/\*/g' *.md
sed -i 's/\$\$\$\$//g' *.md
```

## 常见要求：

- 根据文件_pages中.md的最后修改时间，给文件名添加YYYY-MM-DD才能在博客上显示（rename就行），还要仿照已有的添加frontmatter，tag尽量用和其他类似的，如果有的话。已经有YYYY-MM-DD的就不用了，这种一般frontmatter也都有了。about.md，index这种不要改。
- 在_pages\archive\qq，内容筛选的脚本是现成的，基本不用改。能否根据（大概）2025以来的内容总结出几篇技术记录，放在_pages\Techniques底下？每一篇不太长，就是一篇推送的长度，包括一个或多个话题。尽量按照话题来组织汇总，把比较相关的内容放在一起。内容太多就创建新的md文件。如只看到一个网站但不知道是干啥的，可以访问该网站获取内容。一些文字如果提供了参考网址，要把网址留在下面作为reference，不是文档最下面，就是笔记原位。不要泄露隐私内容。对于这些文档，title: 和第一个# 要加点【笔记整理|2025-07】之类的标签，文件名不用加
