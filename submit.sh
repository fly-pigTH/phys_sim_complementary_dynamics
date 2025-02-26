# Create the submission folder.
mkdir -p submission

# Copy the source code.
cp cpp/material/src/corotated_material.cpp submission/
cp cpp/rig/src/linear_blend_skinning_rig.cpp submission/
cp cpp/sim/src/simulator_basic.cpp submission/
cp cpp/sim/src/simulator_step.cpp submission/

zip -r submission.zip submission