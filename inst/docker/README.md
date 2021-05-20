# ImmuneSignatures2 Docker Image

This directory contains sources for building the docker image used in the ImmuneSpace servers for running all ImmuneSignatures2 code. It can also be used to mimic the server environment if you want to run ImmuneSignatures2 code on your own machine.

## Quickstart

1. Build the image. You must pass in your github personal access token so that the private ImmuneSignatures2 repo can build. The `makefile` has a handy `build` target. It assumes you have `GITHUB_PAT` saved as an environment variable, which is needed to build the ImmuneSignatures2 package, which is in a private repo.

```
gmake build
```

1. Start a container. To run an interactive container:

```
gmake start
```

## How it's set up

Based on a combination of the recommended [LabKey sandboxed docker setup](https://github.com/LabKey/docker-rstudio/tree/develop/images/labkey/rsandbox-ver) and [renv docker setup](https://rstudio.github.io/renv/articles/docker.html)

```
docker
├── Dockerfile
├── README.md
├── install.R
├── makefile
└── renv.lock
```

### `Dockerfile`:

- Starting with [https://github.com/LabKey/docker-rstudio/blob/develop/images/labkey/rsandbox-ver/Dockerfile](https://github.com/LabKey/docker-rstudio/blob/develop/images/labkey/rsandbox-ver/Dockerfile)
  - _add `GITHUB_PAT` environment variable to allow private package installation (TODO: Remove this once repo is public)_
  - install renv ([https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv-1](https://rstudio.github.io/renv/articles/docker.html#creating-docker-images-with-renv-1))
  - copy renv.lock into `/docker/home`
  - install R packages using renv as docker user
    - Currently just sourcing `install.R` which has that as the only line

### `install.R`

```bash
renv::restore(packages = "ImmuneSignatures2", library = "/usr/local/lib/R/site-library")
```

- Install all packages as detailed in `renv.lock` to default site library
  - Install to site library instead of local project library to avoid having to load renv project when running docker container.
- Errors out when you don't do `packages = "ImmuneSignatures2"` because it does not install packages in the right order. TODO: submit issue to renv.
  - This assumes that all required packages are listed in the ImmuneSignatures2 DESCRIPTION file as either imports or suggests.

### `renv.lock`

renv lock file, detailing all package dependencies and versions, and where to download from.

### `makefile`

```makefile
.RECIPEPREFIX = +

start:
+ docker run --rm -ti immunesignatures2

build:
+ docker build -t immunesignatures2:latest --build-arg GITHUB_PAT=${GITHUB_PAT} .
```

- Shortcuts for building image and starting an interactive container.

## Updating `renv.lock` with new dependencies (WIP)

_TODO_

## Installation on ImmuneSpace server

1. [Configure Docker](https://www.notion.so/rglab/Set-up-sandboxed-docker-R-engine-for-ImmSig2-df1e67eeaaff40748983d0d492472ece#9a300ef91df841ea8a5d874142a72752)
2. [Configure labkey](https://www.notion.so/Set-up-Local-ImmuneSpace-753cd9d0df65451e828da9e56f020b2e)
3. `./build.sh` to build the docker image
   1. Make sure you have `GITHUB_PAT` as environment variable
4. Configure labkey to use docker engine [https://www.labkey.org/Documentation/wiki-page.view?name=rsandbox](https://www.labkey.org/Documentation/wiki-page.view?name=rsandbox)
   1. Enable R Docker Sandbox in admin console > experimental features `/admin/experimentalFeatures.view?`
   2. Enable docker engine in admin console > views and scripting `/core/configureReportsAndScripts.view?`
      1. add > new R docker engine
      2. Name: ImmuneSignatures2 R Docker Engine
      3. Docker Image Name: immunesignatures2:4.0.2
      4. Home Directory: /home/docker
      5. "Site Default" Should be unchecked, "Sandboxed" and "Enabled" should be checked.
   3. Configure immunesignatures2 docker engine for IS2 folder [https://www.labkey.org/Documentation/wiki-page.view?name=configureScripting](https://www.labkey.org/Documentation/wiki-page.view?name=configureScripting)
      1. From IS2 folder: Folder > management
      2. Under R config tab:
         1. check "Use folder level R configuration"
         2. Choose "ImmuneSignatures2 R Docker Engine" for reports and pipeline jobs.
         3. Save
      3. If you run a report in IS2 folder, it should use the is2 docker container!
