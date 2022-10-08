# Crescent Moon Visibility Maps Generator

Code is provided to draw crescent visibility maps according to the following criterias:
- [Yallop criteria](https://astro.ukho.gov.uk/download/NAOTN69.pdf),
- [Odeh criteria](https://www.astronomycenter.net/pdf/2006_cri.pdf).

The following features are currently available:
- waxing (evening) crescent visibility bands.
- waning (morning) crescent visibility bands.
- moonset before sunset in red (similarly moonrise after sunrise).

# Output Examples
## Wed, 29 June 2022 (29 ZH 1443H) Evening Crescent
### Yallop
![2022-06-29 (E) Yallop](https://user-images.githubusercontent.com/833473/193407627-e8895f15-7d6f-46c2-9c7f-a770131ad387.png)

For comparison
![2022-06-29 (E) Yallop (HMNAO)](https://user-images.githubusercontent.com/84683703/191850568-3f661abb-74f2-4720-b256-1404d69757cc.jpg)

### Odeh
![2022-06-29 (E) Odeh](https://user-images.githubusercontent.com/833473/193407716-07674584-06c5-47eb-944b-5a6a8ba182bb.png)
  
For comparison
![image](https://user-images.githubusercontent.com/84683703/191850739-bd009136-5e8d-4d0f-ba1d-aac2ace6a564.png)

# Build Instructions
On a Linux machine make sure a C++ compiler (gcc-c++ or g++ package) is available, on macOS install llvm from brew. Make sure imagemagick is installed also, then run `./run.sh`.

On a Windows machine, first install msys2 then run `./build-mingw.sh`

# Credits
- [astronomy-engine](https://github.com/cosinekitty/astronomy/)

# License
MIT
