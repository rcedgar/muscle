#!/bin/bash -e

'c:/Program Files/Microsoft Visual Studio/2022/Community/Msbuild/Current/Bin/MSBuild.exe' muscle.sln /t:muscle /p:Platform=x64 /p:Configuration=Release
