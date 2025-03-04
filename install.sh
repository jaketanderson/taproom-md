#/bin/bash

# Warn the user that their current conda environment will be used
if command -v conda &> /dev/null; then
    # Conda is available
    output=$(conda info | grep "active environment")
    
    # Use cut to split the output using ":" as the delimiter, and keep the text after the colon
    kept_text=$(echo "$output" | cut -d ':' -f 2-)
    
    # Remove leading and trailing whitespace using parameter expansion
    kept_text_trimmed="${kept_text#"${kept_text%%[![:space:]]*}"}"
    kept_text_trimmed="${kept_text_trimmed%"${kept_text_trimmed##*[![:space:]]}"}"
    current_conda_env="$kept_text_trimmed"
else
    # Conda needs to be installed
    echo "The conda command was not available. Please make sure anaconda/miniconda is installed and initialized."
    exit 1
fi

echo "This startup script will install conda packages and taproom in the following conda environment:"
echo "     $current_conda_env"
echo ""
sleep 0.5s

read -p "Is this correct? (y/n): " confirmation
if [ "${confirmation,,}" != "y" ]; then
    exit 1
fi

mamba install -c conda-forge python openmm openff-toolkit openmmforcefields
git clone https://github.com/slochower/host-guest-benchmarks.git

echo "Done, exiting."

exit 0
