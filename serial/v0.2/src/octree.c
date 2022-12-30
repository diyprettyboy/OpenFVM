// Copyright (C) 2004-2005 by OpenFlower Team
// http://sourceforge.net/projects/openflower/
// Licensed under the GNU GPL Version 2.

// Translated from C++ to C by the OpenFVM team in 01/08/2005

#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
#include "octree.h"

void OctCreateNode(oct_node *node, oct_node *parent, int i)
{

	double xmin = parent->xmin;
	double ymin = parent->ymin;
	double zmin = parent->zmin;
	double xmax = parent->xmax;
	double ymax = parent->ymax;
	double zmax = parent->zmax;
	double xmid = parent->xmid;
	double ymid = parent->ymid;
	double zmid = parent->zmid;

	switch(i)
	{
	case Back_Bottom_Left:
		
		node->xmin = xmin;
	    	node->xmax = xmid;
		node->ymin = ymin;
		node->ymax = ymid;
		node->zmin = zmin;
		node->zmax = zmid;
		
		break;

	case Back_Top_Left:
		
		node->xmin = xmin;
		node->xmax = xmid;
		node->ymin = ymin;
		node->ymax = ymid;
		node->zmin = zmid;
		node->zmax = zmax;
		
		break;

	case Front_Bottom_Left:
		
		node->xmin = xmin;
		node->xmax = xmid;
		node->ymin = ymid;
		node->ymax = ymax;
		node->zmin = zmin;
		node->zmax = zmid;
		
		break;
  
	case Front_Top_Left:
		
		node->xmin = xmin;
		node->xmax = xmid;
		node->ymin = ymid;
		node->ymax = ymax;
	    node->zmin = zmid;
		node->zmax = zmax;
		
		break;

	case Back_Bottom_Right:
		
		node->xmin = xmid;
		node->xmax = xmax;
		node->ymin = ymin;
		node->ymax = ymid;
		node->zmin = zmin;
		node->zmax = zmid;
		
		break;

	case Back_Top_Right:
		
		node->xmin = xmid;
		node->xmax = xmax;
		node->ymin = ymin;
		node->ymax = ymid;
		node->zmin = zmid;
		node->zmax = zmax;
		
		break;
  
	case Front_Bottom_Right:
		
		node->xmin = xmid;
		node->xmax = xmax;
		node->ymin = ymid;
		node->ymax = ymax;
		node->zmin = zmin;
		node->zmax = zmid;
		
		break;
  
	case Front_Top_Right:
		
		node->xmin = xmid;
		node->xmax = xmax;
		node->ymin = ymid;
		node->ymax = ymax;
		node->zmin = zmid;
		node->zmax = zmax;
		
		break;
	
	}

	node->xmid = (node->xmin+node->xmax)*0.5;
	node->ymid = (node->ymin+node->ymax)*0.5;
	node->zmid = (node->zmin+node->zmax)*0.5;

}

void OctDestroyNode(oct_node *node)
{
	
	int i;

	for(i = 0; i < 8; i++)
		if(node->nodes[i] != NULL) 
			OctDestroyNode(node->nodes[i]);
		
	free(node);
}


int OctDispatch(oct_node *node, double x, double y, double z)
{

	double xmid = node->xmid;
	double ymid = node->ymid;
	double zmid = node->zmid;

	if( (x<xmid) && (y<ymid) && (z<zmid))
		return Back_Bottom_Left;
  
	if( (x<xmid) && (y<ymid) && !(z<zmid))
		return Back_Top_Left;
  
	if( (x<xmid) && !(y<ymid) && (z<zmid))
		return Front_Bottom_Left;
  
	if( (x<xmid) && !(y<ymid) && !(z<zmid))
		return Front_Top_Left;
  
	if( !(x<xmid) && (y<ymid) && (z<zmid))
		return Back_Bottom_Right;
  
	if( !(x<xmid) && (y<ymid) && !(z<zmid))
		return Back_Top_Right;
  
	if( !(x<xmid) && !(y<ymid) && (z<zmid))
		return Front_Bottom_Right;
  
	if( !(x<xmid) && !(y<ymid) && !(z<zmid))
		return Front_Top_Right;


	printf("\nError: Octree failure\n");
	exit(LOGICAL_ERROR);

}

void OctAddLeaf(oct_node *node)
{
	nbleafs++;
	
	if (nbleafs == 1)
		leafs = malloc(nbleafs * sizeof(oct_node));
	else
		leafs = realloc(leafs, nbleafs * sizeof(oct_node));

	leafs[nbleafs - 1] = *node;
}


void OctCreateRecursive(oct_node *node, oct_data *Tab)
{
	
	int i, j, k, l, m;

	int ok[8];

	int number_of_the_node = 0;

	double xm;
	double xp;
	double ym;
	double yp;
	double zm;
	double zp;

	for (i = 0; i < 8; i++)
	{
		node->nodes[i] = malloc(sizeof(oct_node));
		
		OctCreateNode(node->nodes[i], node, i);
	}
  
	for (i = 0; i < 8; i++) 
	{
		node->nodes[i]->nbentities = 0;
		node->nodes[i]->entities = malloc(sizeof(int));
	}

	for(k = 0; k < node->nbentities; k++)
	{
 
		for (l = 0; l < 8; l++) 
			ok[l] = 0;
		
		i = node->entities[k];
	
		if (i < 0)
		{
			printf("\nError: Octree failure\n");
			exit(LOGICAL_ERROR);
		}

		xm = Tab[i].x-EPSILON;
		xp = Tab[i].x+EPSILON;
		ym = Tab[i].y-EPSILON;
		yp = Tab[i].y+EPSILON;
		zm = Tab[i].z-EPSILON;
		zp = Tab[i].z+EPSILON;
		
		number_of_the_node = OctDispatch(node, xm, ym, zm);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;

		number_of_the_node = OctDispatch(node, xm, ym, zp);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xm, yp, zm);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xm, yp, zp);		
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xp, ym, zm);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xp, ym, zp);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xp, yp, zm);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;
		
		number_of_the_node = OctDispatch(node, xp, yp, zp);
		if(ok[number_of_the_node]==0) ok[number_of_the_node]=1;

		for(m = 0; m < 8; m++)
		{
			if(ok[m] == 1)
			{
				
				node->nodes[m]->nbentities++;
				
				j = node->nodes[m]->nbentities;

				node->nodes[m]->entities = realloc(node->nodes[m]->entities, (j + 1) * sizeof(int));

				node->nodes[m]->entities[j-1] = i;
			}	
		}
	}		

	for(i = 0; i < 8; i++)
	{
  
		if (node->nodes[i]->nbentities > MAXENTITIES)  // Then we need to construct another node starting from this one 
		{
			OctCreateRecursive(node->nodes[i], Tab);
		}
		else
		{
			if (node->nodes[i]->nbentities != 0)  
				OctAddLeaf(node->nodes[i]);
		}
	}

}

void OctCreateOctree(double min[3], double max[3], oct_data *Tab, int nbentities)
{

	int i;

	nbleafs = 0;

	min[0]=min[0]-2.0*EPSILON;
	max[0]=max[0]+2.0*EPSILON;
	min[1]=min[1]-2.0*EPSILON;
	max[1]=max[1]+2.0*EPSILON;
	min[2]=min[2]-2.0*EPSILON;
	max[2]=max[2]+2.0*EPSILON;

	root = malloc(sizeof(oct_node));
	
	root->xmin=min[0];
	root->xmax=max[0];
	root->ymin=min[1];
	root->ymax=max[1];
	root->zmin=min[2];
	root->zmax=max[2];

	root->xmid=(min[0]+max[0])*0.5;
	root->ymid=(min[1]+max[1])*0.5;
	root->zmid=(min[2]+max[2])*0.5;

	root->entities = malloc(nbentities * sizeof(int));
	root->nbentities = nbentities;

	for(i = 0; i < nbentities; i++) 
		root->entities[i] = i;
		
	OctCreateRecursive(root, Tab);

}

void OctDestroyOctree()
{

	if (root != NULL)
	{
		free(root);
		root = NULL;
	}

	if (leafs != NULL)
	{
		free(leafs);
		leafs = NULL;
	}

}

