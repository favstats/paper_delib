# install.packages("webshot")

library(webshot)
install_phantomjs()

webshot(here::here("slides", "slides.html"), "slides/pdf_version.pdf")






# system("/Users/chester/Desktop/infer_workshop/talks/ness-infer/decktape-1.0.0/phantomjs/Users/chester/Desktop/infer_workshop/talks/ness-infer/decktape-1.0.0/decktape.js /Users/chester/Desktop/infer_workshop/talks/ness-infer/slide_deck.html /Users/chester/Desktop/infer_workshop/talks/ness-infer/slide_deck.pdf")

# To R script
# knitr::purl("ness-infer/index.Rmd", "ness-infer/slide_code.R", documentation = 2L)