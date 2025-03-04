# taproom-md

A script to automatically run MD simulations of all complexes in TAPROOM.

*To run:* create a new empty anaconda environment on your system and activate it. Then run `bash install.sh` to install the necessary python libraries and taproom itself. Then all you should have to do is run `python main.py`. If you have issues with CUDA (seeing an error like `CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)`), you may need to change your anaconda environment's version of `cudatoolkit`. For any troubleshooting message me.
