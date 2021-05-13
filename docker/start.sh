# the path to an renv cache on the host machine
RENV_PATHS_CACHE_HOST=/Users/hmiller/Projects/ImmSig2/ImmuneSignatures2/docker/cache

# the path to the cache within the container
RENV_PATHS_CACHE_CONTAINER=/renvcache

# run the container with the host cache mounted in the container
docker run --rm \
    -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
    -e "GITHUB_PAT=${GITHUB_PAT}" \
    -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
    rsandbox:4.0.2 \
    R --vanilla -s -e 'renv::restore()'