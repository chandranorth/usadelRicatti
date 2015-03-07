rem If your compiler is not found, add the path it was installed to in PATH
set PATH=%PATH%;c:\program files\GFortran

python setup.py config_fc --fcompiler=gnu95 --noarch build --compiler=mingw32
python setup.py config install --skip-build
pause

