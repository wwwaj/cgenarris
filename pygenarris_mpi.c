#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stddef.h>
#include "mpi.h"
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "cgenarris_mpi.h"
#include "pygenarris_mpi.h"

#define ZMAX 192
#define VOL_ATTEMPT  100000
#define GRAIN_SIZE 10000

unsigned int *seed;
unsigned int *seed2;
#pragma omp threadprivate(seed)
#pragma omp threadprivate(seed2)
extern float TOL;


//void create_mpi_xtal_type(MPI_Datatype* XTAL_TYPE, int total_atoms);


void mpi_generate_molecular_crystals_with_vdw_cutoff_matrix(
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double tol1, 
	long max_attempts, 
	MPI_Comm world_comm)
{
	float Zp_max=1;
    float volume_mean = volume_mean1;
    float volume_std = volume_std1;
    TOL = tol1;

    if (dim1 != dim2)
	{printf("***ERROR:vdw cutoff matrix is not square***\n"); exit(0);}

    printf("The vdw distance cutoff matrix is :\n");
    for (int i = 0; i < dim1; i++)
    {
    	for (int j= 0; j < dim2; j++)
    	{
    		printf("%f ", *(vdw_matrix + i*dim1 + j) );
    	}
    	printf("\n");
    }

	//Initialise MPI 
	int total_ranks;
    MPI_Comm_size(world_comm, &total_ranks);
    int my_rank;
    MPI_Comm_rank(world_comm, &my_rank);


	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	//variable declarations
	int stop_flag = 0;	// to stop threads if limit is reached
	int success_flag = 0;
	int counter = 0;	//counts number of structures
	int spg_index = 0;	//space group to be generated
	FILE *out_file;		//file to output geometries
	if (my_rank == 0)
	{
		out_file = fopen("geometry.out","w");
		if(!out_file)		//check permissions
		{
			printf("***ERROR: cannot create geometry.out \n");
			exit(0);
		}
		//fprintf(out_file, "my_rank=%d\n", my_rank);
	}


	//random number seeding, different seeds for different threads
	seed = (unsigned int*)malloc(sizeof(unsigned int)); //seed for uniform gen
	seed2 = (unsigned int*)malloc(sizeof(unsigned int)); //seed for random
	*seed += my_rank*7 + seed_shift*13; //some random seed private for each threads
	*seed2 = my_rank*17 + seed_shift*11;
	init_genrand(abs(*seed));
	
	//storing information for compatible space groups
	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;	//storing number molecule axes
	float *mol_axes; //storing possible molecular axes
	
	//variable declarartion	
	molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));//dummy crystal
	//float volume_std;	//standard dev for volumes
	//float volume_mean;	//mean volume
	float sr = -1;			//specific radius proportion for structure checks
						//see paper for definition
	//float Zp_max;		//Z'' . not implemented
	float volume;		//random volume used of generation
	//int Z;				//multiplicity of general position
	//int num_structures;	//num of structures
	int spg;			//space group attempted
	//long max_attempts;	//max attempts per space group

	//read input from file, read molecular structure from geometry,Z, Zp
	if(my_rank == 0)
	{
		int len = 100;
		char name[len];
		gethostname(name, len);
		printf("PARALLELIZATION INFO:\n");
		printf("---------------------------\n");
		printf("Using MPI for parallelization.\n");
		printf("cgenarris is running on %d processes. \n", total_ranks);
	}

	MPI_Barrier(world_comm);
	for(int i = 0; i < total_ranks; i++)
	{	
		if(i == my_rank)
			printf("Rank number %d out of %d reporting.\n", my_rank+1, total_ranks);
    	MPI_Barrier(world_comm);
    }

	if(my_rank == 0)
	{	
		printf("---------------------------\n");
		printf("\n");
	}


    read_geometry(mol);				//read molecule from geometry.in
    /*
    read_control(&num_structures,
				 &Z,
				 &Zp_max,
				 &volume_mean, 
				 &volume_std,
				 &sr,
				 &max_attempts);	//get settings
	*/
	
	//recenter molecule to origin
	recenter_molecule(mol);

	if (my_rank == 0)
	{
		print_input_geometry(mol);
		print_input_settings(&num_structures,
							 &Z,
							 &Zp_max,
							 &volume_mean, 
							 &volume_std,
							 &sr,
							 &max_attempts);
	}

	MPI_Barrier(world_comm);

    //inititalise volume
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N); //allcate memory
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	int total_atoms = mol->num_of_atoms * Z;

	//create xtal type structure for MPI
	//MPI_Datatype XTAL_TYPE;
	//create_mpi_xtal_type(&XTAL_TYPE, mol->num_of_atoms*Z);
	//MPI_Type_commit(&XTAL_TYPE);

	
	if(my_rank == 0)
	{
		printf("COMPATIBLE SPACE GROUP INFO:\n");
		printf("-----------------------------------------------------\n");
		printf("Detecting compatible space groups"
			   " from molecule symmetries...\n");
	}

	//compatible space groups are stored in compatible_spg. see the 
	//object definition in spg_generation.h
	//every thread has its own copy
	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  &mol_axes,
								  &num_axes,
								  my_rank+1);
			
	if(my_rank == 0)
	{
		printf("\n");
		printf("Total compatible space groups = %d\n", num_compatible_spg);
		printf("Number of molecular axes = %d \n", num_axes);
		printf("-----------------------------------------------------\n\n");
		printf("Starting generation...\n\n");
		sleep(1);
	}
	
	MPI_Barrier(world_comm); // wait for other friends to join
	
	/*deprecated
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	*/
	
	
	while( spg_index < num_compatible_spg )
	{
		MPI_Barrier(world_comm);
		//counter counts the number of structures generated
		spg = compatible_spg[spg_index].spg; //pick a spg 

		//print attempted space group 
		if (my_rank == 0)
			{printf("Attempting to generate spacegroup number %d....\n", spg);}

		while( counter < num_structures )
		{
			int verdict = 0; //for structure check
			int i = 0; 		 //counts attempts for spg
			//attempts for an spg.
			stop_flag = 0; 
			for(; i < max_attempts/total_ranks; i = i + GRAIN_SIZE) 
			{
				int j = 0;
				success_flag = 0;
				for(; j < GRAIN_SIZE; j++)
				{
					//generate
					int result = generate_crystal(random_crystal,
									 mol,
									 volume,
									 Z,
									 Zp_max,
									 spg,
									 compatible_spg,
									 num_compatible_spg,
									 spg_index,
									 mol_axes,
									 num_axes);
					
					//alignment failure
					if(!result)
						continue;
					
					//check if molecules are too close with sr	    
					verdict = check_structure_with_vdw_matrix(*random_crystal, vdw_matrix, dim1, dim2);    
					
					//reset volume after volume attempts
					if( (i+j) % VOL_ATTEMPT == 0 && i+j != 0)
					{
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
						if(my_rank == 0)
							printf("#Rank 0: Completed %d attempts.\n", (i+j)*total_ranks);
						fflush(stdout);
					}

										//if generation is successful
					if (verdict == 1 )
					{
						random_crystal->spg = spg;
						printf("#Rank %d: Generation successful.\n",
							my_rank);					
						printf("#attempted space group = %d,\n"
								"#unit cell volume (cubic Angstrom)= %f.\n",  
								spg,
								volume);
						int spglib_spg = detect_spg_using_spglib(random_crystal);
						printf("#SPGLIB detected space group = %d\n\n",
											                     spglib_spg);
						break;
						//verdict = 1;
					}


				}//end of GRAIN loop

				int found_poll[total_ranks];
				//printf("verdict = %d\n", verdict);
				MPI_Gather(&verdict, 1, MPI_INT, &found_poll, 1, MPI_INT, 0, world_comm);
				if (my_rank == 0)
				{
					//print the structure generated by root first to outfile
					if(verdict)
					{
						if (counter < num_structures)
						{							
							print_crystal2file(random_crystal, out_file);
							counter++;
							success_flag = 1;
						}
						else
							stop_flag = 1;
						
					}
					//get and print structures from other ranks
					for(int rank = 1; rank < total_ranks; rank++)
					{
						if (found_poll[rank] == 1)
						{
							receive_xtal(world_comm, rank, random_crystal, total_atoms);
							//print_crystal(random_crystal);
							success_flag = 1; 
							if(counter < num_structures)
							{	
								print_crystal2file(random_crystal, out_file);
								printf("#Rank %d: Generation successful.\n", rank);
								counter++;
							}
							else
								stop_flag = 1;
						}
					}
				}
				//all other ranks send the crystal to root rank
				else
				{
					if (verdict == 1)
					{
						send_xtal(world_comm, 0, random_crystal, total_atoms);
						i = 0;
					}	
				}

				MPI_Bcast(&success_flag, 1, MPI_INT, 0, world_comm);
				MPI_Bcast(&stop_flag, 1, MPI_INT, 0, world_comm);

				if (success_flag)
				{
					i = 0;
					do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				}
				if (stop_flag)
					break;


			}//end of attempt loop	
			
			//if max limit is reached if some rank hit the limit
			if (i >= max_attempts/total_ranks)
			{	
				do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				if (my_rank== 0)
				{	
					printf("**WARNING: generation failed for space group = %d "
							"after %ld attempts. \n",
							spg,
							max_attempts);		
					fflush(stdout);
					//print_crystal(random_crystal);
				}
				counter = num_structures + 1;
				MPI_Barrier(world_comm);
			}

			if(stop_flag)
				break;	
			
		}//end of numof structures whileloop
	
		//move to next spacegroup
		counter = 0;
		spg_index++;
		if (my_rank == 0)
		{
			 printf("#space group counter reset. Moving to next space group...\n\n");
			 fflush(stdout);
		} 
		MPI_Barrier(world_comm);
		
	}//end of spg while loop

	if(my_rank == 0)
		fclose(out_file);

	if(my_rank == 0)
		printf("Generation completed.\nHave a nice day!!\n");

}