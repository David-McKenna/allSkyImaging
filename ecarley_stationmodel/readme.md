# Known Issues
* Matplotlib hates giving you control of the z-buffer. Projections will often cover the array scatter points or sun itself.
* Some projections will exit the axis limits; extending the axis limits to accomate this gives uglier results and would require rescaling all the axis each time, not worth the hassle.