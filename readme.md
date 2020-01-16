# canopylife

Code repository to run analyses and generate figures and manuscript for Nitta *et al.* "Life in the canopy: Community trait assessments reveal substantial functional diversity among fern epiphytes". FIXME: ADD COMPLETE CITATION WHEN AVAILABLE

All code is in [R](https://cran.r-project.org/). The [drake package](https://ropensci.github.io/drake/) is used to manage the workflow. To run all analyses and generate the manuscript, [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository), download the [data](#data), and run `make.R`.

## Data

All data are stored on [Dryad](https://datadryad.org/). After cloning the repository, download the following sets of zipped data to the `data/` folder (click on the "Download Dataset" button for each link below):

- https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g

- FIXME: ADD DRYAD LINK WHEN AVAILABLE

The code needs the zipped data files to run, so if you unzip them yourself to inspect them, be sure to keep the original, zipped versions in `data/`.

## Reproducible analysis with Docker

`make.R` requires various software packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/canopylife) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/).

Clone this repository and download the data as described above.

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

You will see the targets being built by `drake`, and the final manuscript should be compiled at the end as `manuscript.pdf` (for easy viewing) and `manuscript.docx` (for journal submission) in the `results/` folder. Other results files (image files, SI, etc.) will also be output to `results/`.

For submission, a very small number of manual tweaks to `manuscript.docx` were made where indicated with `FIXME` comments in `ms/manuscript.Rmd`.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```
