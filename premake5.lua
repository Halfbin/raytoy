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
    symbols "On"

  filter "configurations:release"
    defines { "NDEBUG" }
    optimize "On"

  filter "configurations:profile"
    defines { "NDEBUG" }
    optimize "On"
    symbols "On"
    buildoptions { "-pg" }
    linkoptions { "-pg" }

  filter "action:gmake" -- FIXME: this should be toolset:gcc, but toolset: is broken in premake5 as of 2015-09-01
  --filter "toolset:gcc"
    buildoptions { "-std=c++14" }

