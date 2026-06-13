import struct

def inspect_bin(file_path):
    print(f"--- Inspecting Binary File: {file_path} ---")
    try:
        with open(file_path, "rb") as f:
            data = f.read()
            
        # Scan the binary chunk for printable ASCII variable names
        # like 'X', 'Y', 'F_1', 'X_1' etc.
        import re
        strings = re.findall(b'[A-Za-z0-9_]{1,20}', data)
        
        unique_params = sorted(list(set([s.decode('ascii') for s in strings])))
        
        print("\nFound the following Parameter symbols inside the file:")
        for param in unique_params:
            if any(c.isalpha() for c in param): # filter out pure numbers
                print(f"  -> {param}")
                
    except FileNotFoundError:
        print(f"Error: Could not find file at {file_path}")

if __name__ == "__main__":
    inspect_bin("./build/generator/derivatives.bin") # Change path if it's inside your build/ folder