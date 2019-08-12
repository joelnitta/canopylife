# canopylife

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Constrained, yet diverse: Phylogenetic and morphological filtering in epiphytic ferns".

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow. To run all analyses and generate the manuscript, simply [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `make.R`.

## Reproducible analysis with Docker

`make.R` requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/canopylife) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/) and clone this repository.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Enter the container:

```
docker exec -it canopylife_analysis_1 bash
```

Inside the container, run `make.R`:

```
Rscript make.R
```

You will see the targets being built by `drake`, and the final manuscript should be compiled at the end as `manuscript.pdf` in the `ms` folder. Other figure pdfs, tables (in rich-text format), and the SI (.doc format) will also be compiled.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```
