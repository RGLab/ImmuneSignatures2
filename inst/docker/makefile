.RECIPEPREFIX = +

GITHUB_PAT=do_not_commit_your_github_pat

start:
+ docker run --rm -ti immunesignatures2

build: 
+ docker build -t immunesignatures2:latest --build-arg GITHUB_PAT=${GITHUB_PAT} .
