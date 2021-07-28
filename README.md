# MultiTileBlueNoise
Generates blue noise with an expansion on https://github.com/Atrix256/VoidAndCluster ; also included a library of generated blue noise textures.

You are probably going to want the blue noise textures rather than the generator, to be honest. The generator is
not much different from the VoidAndCluster repo it is based on. The generated textures are... well, more special.
They are in the `results` folder, and include "omni-tiling" blue noise textures, which are like Wang tiles but
have no constraints for what other textures they may be placed next to (they tile with any texture in their
group). There are also triangular-mapped versions of these. The best omni-tiling textures are in the `omni_equal`
and `omni_tri_equal` directories; these are equivalent to normal blue noise in how they distribute color/byte
values, and have 32x32, 64x64, and 128x128 textures. The `omni` and `omni_tri` directories are also usable, but
aren't as even with how frequently colors appear. The `standard` and `tri` folders have normal blue noise
textures that don't have the omni-tiling property, though they do have the correct color frequencies. The folders
in `results` have small files with names that say what the sigma is for that directory; it's always 1.9 now, but
the older blue noise textures in `results` (not a subdirectory of it) have a higher sigma that I don't remember.

All of the textures in `results` are made with the BlueNoise classes in
[SquidSquad's tests](https://github.com/yellowstonegames/SquidSquad/tree/main/squidworld/src/test/java/com/github/yellowstonegames/world);
that's Java code. It's short but quite complex. The omni-tiling generators work by using a modified version of Ulichney's
void-and-cluster technique, based off of
[Bart Wronski's beautifully simple generator](https://bartwronski.com/2021/04/21/superfast-void-and-cluster-blue-noise-in-python-numpy-jax/).
They split the parent blue noise into 64 or 16 sectors (each will be used to make a different blue noise texture), and where void-and-cluster
(in Wronski's version) normally propagates "energy" into the surrounding area, wrapping at the parent's edges, when this propagates energy
across any sector edge, it wraps across every corresponding opposite edge of a sector. The omni_equal varieties also track how many times
each sector has been assigned a value, and will avoid placing more energy into a sector that has already been given the full allotment.
