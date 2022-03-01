library(bookdown)
library(dplyr)
library(kableExtra)

# Bookdown resources https://bookdown.org/yihui/bookdown/github.html

# Confirms the nojekyll is created
file.create('\\\\142.244.87.176/abmisc/github-repos-data/LandFacets/docs/.nojekyll')

# Render bookdown
bookdown::render_book(input = "bookdown/", 
                      output_format = "bookdown::gitbook",
                      output_dir = "\\\\142.244.87.176/abmisc/github-repos-data/LandFacets/docs/")




