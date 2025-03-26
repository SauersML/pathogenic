#!/bin/bash
set -e

# Print banner
echo "====================================="
echo "Pathogenic Variant Finder Installer"
echo "====================================="
echo

# Check if Rust is installed
if ! command -v cargo &> /dev/null; then
    echo "Error: Rust and Cargo are required but not installed."
    echo "Please install Rust from https://rustup.rs/ and try again."
    exit 1
fi

echo "Building release version..."
cargo build --release

BINARY_PATH="$(pwd)/target/release/pathogenic_variant_finder"
echo "Binary built at: $BINARY_PATH"

# Ask user for installation preference
echo
echo "Installation options:"
echo "1) Create symlink in /usr/local/bin (requires sudo, recommended)"
echo "2) Create symlink in ~/.local/bin (no sudo required)"
echo "3) Skip installation (just build)"
echo
read -p "Select an option [1-3]: " choice

case $choice in
    1)
        if [ ! -d "/usr/local/bin" ]; then
            sudo mkdir -p /usr/local/bin
        fi
        echo "Creating symlink in /usr/local/bin..."
        sudo ln -sf "$BINARY_PATH" /usr/local/bin/pathogenic
        echo "Symlink created: /usr/local/bin/pathogenic"
        echo "You can now run the tool using the 'pathogenic' command."
        ;;
    2)
        if [ ! -d "$HOME/.local/bin" ]; then
            mkdir -p "$HOME/.local/bin"
        fi
        echo "Creating symlink in ~/.local/bin..."
        ln -sf "$BINARY_PATH" "$HOME/.local/bin/pathogenic"
        echo "Symlink created: $HOME/.local/bin/pathogenic"
        
        # Check if ~/.local/bin is in PATH
        if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
            echo "Note: $HOME/.local/bin is not in your PATH."
            echo "Add the following line to your ~/.bashrc or ~/.zshrc:"
            echo "  export PATH=\"\$HOME/.local/bin:\$PATH\""
            echo "Then restart your terminal or run 'source ~/.bashrc' (or ~/.zshrc)"
        else
            echo "You can now run the tool using the 'pathogenic' command."
        fi
        ;;
    3)
        echo "Installation skipped. To run the tool, use:"
        echo "  $BINARY_PATH [arguments]"
        ;;
    *)
        echo "Invalid option. Installation skipped."
        echo "You can manually run the tool with:"
        echo "  $BINARY_PATH [arguments]"
        ;;
esac

echo
echo "Installation complete!"
echo "Run 'pathogenic --help' for usage information." 