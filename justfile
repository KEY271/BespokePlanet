root := justfile_directory()

# Extra compiler flags on top of fpm's release profile: automatic arrays on the
# stack and the instruction set of the build machine.
fflags := "-fstack-arrays -march=native"

# Run the complete Fortran test suite in release mode.
test:
    cd "{{ root }}/core" && fpm test --profile release --flag "{{ fflags }}" --link-flag "$(pkg-config --libs-only-L fftw3)"

# Run the model in release mode (equation: shallow-water, barotropic, dry, held-suarez,
# radiation, slab-ocean, moist, or all).
run equation="dry":
    cd "{{ root }}/core" && fpm run --profile release --flag "{{ fflags }}" --link-flag "$(pkg-config --libs-only-L fftw3)" -- "{{ equation }}"

# Start the visualizer (defaults to http://127.0.0.1:8000).
visualizer port="8000":
    cd "{{ root }}/viz" && if [ ! -f node_modules/three/build/three.module.js ]; then npm install; fi
    cd "{{ root }}/viz" && python3 server.py --output "{{ root }}/output" --port "{{ port }}"

alias viz := visualizer
