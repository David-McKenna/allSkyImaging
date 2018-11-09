uvBeamOnly -- UNFINISHED, FFT issues. Rest should be fine.

# Known Issues
* If run outside of an IPython session, the 3D plotting will be unglodly slow. Use the run -i ./scripy.py <flags> syntax to call it.
* Matplotlib hates giving you control of the z-buffer. Projections will often cover the array scatter points or sun itself.
* Some projections will exit the axis limits; extending the axis limits to accomate this gives uglier results and would require rescaling all the axis each time, not worth the hassle.