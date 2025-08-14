#!/usr/bin/env python3
"""
Test the specific import that's failing.
"""

import sys
from pathlib import Path

def test_import_step_by_step():
    """Test each step of the import process."""
    print("Testing EasyMD imports step by step...")
    print("=" * 50)
    
    # Add src to path
    src_dir = Path.cwd() / "src"
    if src_dir.exists():
        sys.path.insert(0, str(src_dir))
        print(f"✅ Added {src_dir} to Python path")
    else:
        print(f"❌ {src_dir} does not exist")
        return False
    
    # Test step 1: Import EasyMD package
    try:
        import EasyMD
        print("✅ Step 1: import EasyMD - SUCCESS")
        print(f"   EasyMD location: {EasyMD.__file__}")
    except Exception as e:
        print(f"❌ Step 1: import EasyMD - FAILED: {e}")
        return False
    
    # Test step 2: Import sysGenerator package
    try:
        import EasyMD.sysGenerator
        print("✅ Step 2: import EasyMD.sysGenerator - SUCCESS")
        print(f"   sysGenerator location: {EasyMD.sysGenerator.__file__}")
    except Exception as e:
        print(f"❌ Step 2: import EasyMD.sysGenerator - FAILED: {e}")
        return False
    
    # Test step 3: Import sysGenerator module
    try:
        import EasyMD.sysGenerator.sysGenerator
        print("✅ Step 3: import EasyMD.sysGenerator.sysGenerator - SUCCESS")
        print(f"   sysGenerator.sysGenerator location: {EasyMD.sysGenerator.sysGenerator.__file__}")
    except Exception as e:
        print(f"❌ Step 3: import EasyMD.sysGenerator.sysGenerator - FAILED: {e}")
        return False
    
    # Test step 4: Import the class
    try:
        from EasyMD.sysGenerator.sysGenerator import SysGenerator
        print("✅ Step 4: from EasyMD.sysGenerator.sysGenerator import SysGenerator - SUCCESS")
        print(f"   SysGenerator class: {SysGenerator}")
    except Exception as e:
        print(f"❌ Step 4: from EasyMD.sysGenerator.sysGenerator import SysGenerator - FAILED: {e}")
        return False
    
    # Test step 5: Import via package __init__
    try:
        from EasyMD.sysGenerator import SysGenerator as SysGen2
        print("✅ Step 5: from EasyMD.sysGenerator import SysGenerator - SUCCESS")
        print(f"   SysGenerator via __init__: {SysGen2}")
    except Exception as e:
        print(f"❌ Step 5: from EasyMD.sysGenerator import SysGenerator - FAILED: {e}")
        return False
    
    print("\n🎉 All import steps successful!")
    return True

def check_file_contents():
    """Check the contents of key files."""
    print("\nChecking file contents...")
    print("=" * 30)
    
    files_to_check = [
        "src/EasyMD/__init__.py",
        "src/EasyMD/sysGenerator/__init__.py",
    ]
    
    for file_path in files_to_check:
        path = Path(file_path)
        if path.exists():
            print(f"\n📄 {file_path}:")
            try:
                with open(path, 'r') as f:
                    content = f.read().strip()
                    if content:
                        print(f"   {content}")
                    else:
                        print("   (empty file)")
            except Exception as e:
                print(f"   Error reading file: {e}")
        else:
            print(f"\n❌ {file_path}: File not found")

if __name__ == "__main__":
    print("EasyMD Specific Import Test")
    print("This tests the exact import that's failing")
    print()
    
    check_file_contents()
    success = test_import_step_by_step()
    
    if not success:
        print("\n🔧 Try these solutions:")
        print("1. Run: python fix_installation.py")
        print("2. Run: python diagnose_installation.py")
        print("3. Use: python -m EasyMD instead of easymd command")