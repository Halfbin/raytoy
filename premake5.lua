workspace "raytoy"
  configurations { "debug", "release", "profile" }

project "raytoy"
  kind "ConsoleApp"
  language "C++"
  targetdir "%{cfg.buildcfg}/bin"

  files { "src/**.hpp", "src/**.cpp" }

  -- includedirs { "rk-core/include", "rk-math/include" }
  -- links { "c++abi", "c++" }
  links { "pthread", "SDL2" }

  warnings "Extra"

  filter "configurations:debug"
    defines { "DEBUG" }
  --buildoptions { "-Og" }
    symbols "On"

  filter "configurations:release"
    defines { "NDEBUG" }
  --optimize "Speed"
    buildoptions { "-flto", "-O3", "-march=core-avx2" }

  filter "configurations:profile"
    defines { "NDEBUG" }
    symbols "On"
    buildoptions { "-flto", "-O3", "-pg", "-no-pie" }
    linkoptions { "-flto", "-O3", "-pg", "-no-pie" }

  buildoptions { "-std=c++14" }

