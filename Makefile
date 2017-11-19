md:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

site:
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

cov:
	Rscrip -e "covr::codecov()"

