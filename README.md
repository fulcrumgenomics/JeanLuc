[![Build Status](https://travis-ci.org/fulcrumgenomics/JeanLuc.svg?branch=master)](https://travis-ci.org/fulcrumgenomics/JeanLuc)


# This project is no longer in active development.

Please see [fgbio](https://github.com/fulcrumgenomics/fgbio) which replaces Jean Luc.

# JeanLuc

## Java Requirements
Java SE 8 is required.

## Building JeanLuc

JeanLuc uses gradle to build. There are two steps involved:

1. Install "gradle" (ideally using `brew install gradle` with [Homebrew](http://brew.sh/), or from http://gradle.org/gradle-download/)
2. In this directory (jeanluc) run `gradle assemble`

You can also
* Run a clean build with `gradle clean assemble`
* Run unit tests with `gradle [clean] test` (you can view a test report at `./build/reports/tests/index.html`)
* Run a single unit test class with `gradle -Dtest.single=[pattern] test` where pattern is a class name or partial class name

After running a build you can find the resulting JAR file in `./build/libs/jeanluc-1.0.0.jar`

## Configuring IntellIJ for JeanLuc

IntelliJ natively understands gradle builds, which is very helpful because it can import the project automatically, and also build it using the gradle build file.  To setup IntelliJ do the following:

1. Open IntelliJ
2. From the menus File -> New -> Project From Existing Sources
3. Use the file navigator to find and select `build.gradle` in this directory and hit OK
4. Select "Use local gradle distribution" and give it the path to your gradle home (e.g. `/usr/local/Cellar/gradle/2.10/libexec`)
5. Ensure the "Project format" is set to `.idea` and click OK
6. Done!
