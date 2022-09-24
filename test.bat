@echo off

set diff="C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\Common7\IDE\CommonExtensions\Microsoft\TeamFoundation\Team Explorer\Git\usr\bin\diff.exe"

out\build\x86-Debug\pipefoil.exe test 1 | %diff% - testfiles\test1.dat
out\build\x86-Debug\pipefoil.exe test 2 | %diff% - testfiles\test2.dat
out\build\x86-Debug\pipefoil.exe test 3 | %diff% - testfiles\test3.dat
out\build\x86-Debug\pipefoil.exe test 4 | %diff% - testfiles\test4.dat
out\build\x86-Debug\pipefoil.exe test 5 | %diff% - testfiles\test5.dat
out\build\x86-Debug\pipefoil.exe 160 4.7 50 1.5 20 | %diff% - testfiles\test5.dat
