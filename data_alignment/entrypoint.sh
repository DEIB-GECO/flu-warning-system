#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# This script is the container's ENTRYPOINT.
# It receives all arguments passed to `docker-compose run`.

# Example: If the first argument is "shell", open a bash shell.
# Otherwise, treat the arguments as a command to execute.

if [ "$1" = "H5N1" ]; then
    # Shift to remove the "shell" argument, so remaining args don't interfere with bash
    shift
    exec /app/H5N1.sh
elif [ "$1" = "H1N1" ]; then
    shift # Remove "script1" from arguments
    exec /app/H1N1.sh
else
    # Default behavior: treat the arguments as a command to execute directly
    echo "Execution aborted because of the mandatory argument 'H1N1' or 'H5N1' is missing."
fi