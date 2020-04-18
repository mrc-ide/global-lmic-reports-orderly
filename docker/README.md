The build is broken into a few steps in order to exploit caching and scheduled builds.

There's a fat base image that does not need building that often.  We'll do this on a trigger from the repository, but not on a scheduled build
