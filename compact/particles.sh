#!/bin/bash

# define particles here 
particle[1]="mu-"
particle[2]="e-"
particle[3]="gamma"
particle[4]="pi-"
particle[5]="pi0"
particle[6]="neutron"
particle[7]="proton"
particle[8]="kaon-"

Nparticles=${#particle[@]}
echo "Nr of particles = $Nparticles"
