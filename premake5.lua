workspace "raytoy"
  configurations { "debug", "release", "profile" }

project "raytoy"
  kind "ConsoleApp"
  language "C++"
  targetdir "%{cfg.buildcfg}/bin"

  files { "src/**.hpp", "src/**.cpp" }

  -- includedirs { "rk-core/include", "rk-math/include" }
  -- links { "c++abi", "c++" }
  links { "pthread" }

  warnings "Extra"

  filter "configurations:debug"
    defines { "DEBUG" }
    buildoptions { "-Og" }
    symbols "On"

  filter "configurations:release"
    defines { "NDEBUG" }
    --optimize "Speed"
    buildoptions { "-flto", "-Ofast" }

  filter "configurations:profile"
    defines { "NDEBUG" }
    optimize "On"
    symbols "On"
    buildoptions { "-pg" }
    linkoptions { "-pg" }

  buildoptions { "-std=c++14" }

