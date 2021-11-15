.RECIPEPREFIX = +
# .RECIPEPREFIX allows you to use something other than tabs.
vignette_path = $(realpath ../../vignettes)
data_raw_path = $(realpath ../../data-raw)

start:
+ docker run --rm -ti \
  -v ${vignette_path}:/home/docker/vignettes \
  -v ${data_raw_path}:/home/docker/data-raw \
  immunesignatures2

build:
+ docker build -t immunesignatures2:latest .
