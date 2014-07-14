/*
 *
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

//Particles
#define Si 0
#define Xe 1
#define XeSi 2
#define C 3
#define H 4
#define SiMass 28.085
#define XeMass 131.293

#define RAND_MAX 0x7fff

//Parameters for Si-Si and Si-Si-Si
#define A_Si 7.049556277
#define B_Si 0.6022245584
#define lamda_Si 21.0
#define gama_Si 1.2
#define sigma_Si 2.0951
#define space_Si 5.4309497784624312826814139567494 //(pow(2.0,(1.0/6.0))*sigma_Si*(4/sqrt(3)))
#define epsilon_Si 0.02092
#define a_Si 1.8
#define si_Cluster 22.0f

//Parameters for Xe-Xe
#define sigma_Xe_Xe 4.05
#define epsilon_Xe_Xe 0.00021
#define space_Xe 6.4289742604712078727444078390528 //2^(1/6)*sigma_Xe*(???)
#define xe_Cluster 22.0f

//parameters for Si-Xe and Xe-Si
#define sigma_Si_Xe 3.9
#define epsilon_Si_Xe 0.00016
#define si_xe_Cluster 8.0
#define a_Si_Xe 2000.0f

// Device Constants
#define GEAR1			0.15833333333333333333333333333333
#define GEAR2			0.75
#define GEAR3			0.5
#define GEAR4			0.08333333333333333333333333333333

#define rb 8.314e-7; // used for temperature

#endif	//CONSTANTS_H_