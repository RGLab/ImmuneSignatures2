# Don't use package cache. https://rstudio.github.io/renv/articles/renv.html#cache-1
# Unneccesary since you only install packages once, when building container. 
# Default cache adds lots of extra files under home directory which makes it really slow when 
# running in LabKey. 
renv::settings$use.cache(FALSE)

# Install packages: renv does not resolve dependencies correctly if you don't do 
# packages = "ImmuneSignatures2". 
# This assumes that ImmuneSignatures2 package lists all relevant additional dependencies in 
# DESCRIPTION file. 

# For simplicity, install packages into default site library, so no need to worry about 
# juggling libPaths. Not recommended outside of containerized R installations! 
renv::restore(packages = "ImmuneSignatures2", library = "/usr/local/lib/R/site-library") 
