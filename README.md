# MultiTileBlueNoise
Generates blue noise with an expansion on https://github.com/Atrix256/VoidAndCluster ; also included a library of generated blue noise textures.

You are probably going to want the blue noise textures rather than the generator, to be honest. The generator
is not especially fast and doesn't do anything special. The generated textures are... well, more special. They include
"omni-tiling" blue noise textures, which are like Wang tiles but have no constraints for what other textures they
may be placed next to. There are also triangular-mapped versions of these. Optimally, you would have the textures
you want in the `omni_equal` or `omni_tri_equal` directories; these are equivalent to normal blue noise in how they
distribute color/byte values, but these folders aren't complete yet, and only have 32x32 textures. The `omni` and
`omni_tri` directories are more complete currently, but aren't as even with how frequently colors appear.
`standard` and `tri` have normal blue noise textures that don't have the omni-tiling property, though they do
have the correct color frequencies. All of these textures are made with the BlueNoise classes in
[SquidSquad's tests](https://github.com/yellowstonegames/SquidSquad/tree/main/squidworld/src/test/java/com/github/yellowstonegames/world);
that's Java code. It's short but somewhat complex. 
