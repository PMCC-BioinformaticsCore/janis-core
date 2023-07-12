


"""
This exists so that we can do tests on tools, types, workflows which are in the  janis-bioinformatics & janis-unix packages.

This is not current with janis-bioinformatics / janis-unix. 

janis-bioinformatics & janis-unix have their own tests - the ones here are 
just so that we can test janis-core functionality with things that at one time were in 
the above packages. 

Basically just a way to not need janis-bioinformatics & janis-unix as dependencies,
but can still test legacy tools, types, workflows. 

Previously, this was the situation
- janis-core imports janis-bioinformatics
- janis-core imports janis-unix
- janis-unix imports janis-core
- janis-bioinformatics imports janis-core
- janis-bioinformatics imports janis-unix. 

We are trying to decouple these seperate packages. 

"""