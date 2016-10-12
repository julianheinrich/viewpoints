# Periodic Orbit Example
This example demonstrates how to create a camera orbit around a fixed axis. The resulting orbit will pass both low and high viewpoint-entropy views for the given selection. This ensures that both the composition of complexes (which often correspond to low-entropy viewpoints) as well as the overall structure (high entropy viewpoints) should be visible during the rotation:

<img src='orbit.gif'></img>

This example was created using the following commands:

```python
fetch 2d1s
as cartoon
colorByAquaria
as spheres, organic
orbit all, by=ss, n=100, width=512, height=512
```