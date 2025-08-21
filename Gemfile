# Gemfile (Corrected Version)

source "https://rubygems.org"

# Specifies the Jekyll version for your project.
gem "jekyll", "~> 4.3.3"

# Switch to the default, stable Jekyll theme to avoid Sass issues.
gem "minima", "~> 2.5"

# Specifies plugins for your site.
group :jekyll_plugins do
  gem "jekyll-feed", "~> 0.12"
end

# Windows-specific gems for timezone data and performance.
platforms :mingw, :x64_mingw, :mswin, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
  gem "wdm", "~> 0.1.1"
end

# JRuby-specific gem lock.
gem "http_parser.rb", "~> 0.6.0", :platforms => [:jruby]