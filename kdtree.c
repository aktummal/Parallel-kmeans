#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

int cluster_counter = 0;

void kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid);
void bipartition(int dim, int i0, int im, double *data, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid);
void search_kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double *query, double *result_pt);

void main(){
	int dim=0, ndata=0, k=0, i=0, j=0;
	printf("Enter the number of dimensions: ");
	scanf("%d",&dim); // Dimensions input
	printf("Enter the number of points: ");
	scanf("%d",&ndata); //No. of points input
	printf("Enter the number of clusters: ");
	scanf("%d",&k); //No. of clusters input
	double *data; //Data array to store the points
	data = (double*) malloc((dim*ndata) * sizeof(double));
	int *cluster_size; //To store the sizes of the clusters
	cluster_size = (int*) malloc((2*k-1) * sizeof(int));
	int *cluster_start; //To store the starting indexes of the clusters
	cluster_start = (int*) malloc((2*k-1) * sizeof(int));
	double **cluster_bdry; //To store the min and max values of each dimension for each cluster
	cluster_bdry = (double**) malloc((2*k-1) * sizeof(double*));
	for(i=0;i<(2*k-1);i++)
		cluster_bdry[i] = (double*) malloc((dim*2) * sizeof(double));
	double **cluster_centroid; //To store the centroid of each cluster
	cluster_centroid = (double**) malloc((2*k-1) * sizeof(double*));
	for(i=0;i<(2*k-1);i++)
		cluster_centroid[i] = (double*) malloc(dim * sizeof(double));
	for(i=0;i<2*k-1;i++){ //Initializing the values to -1 initially
		cluster_size[i]=-1;
		cluster_start[i]=-1;
	}
	for(i=0;i<2*k-1;i++) //Initializing the values to 0
		for(j=0;j<dim*2;j++)
			cluster_bdry[i][j]=0.0;
	for(i=0;i<2*k-1;i++) //Initializing the values to 0
		for(j=0;j<dim;j++)
			cluster_centroid[i][j]=0.0;
	kdtree(dim,ndata,data,k,cluster_size,cluster_start,cluster_bdry,cluster_centroid); //Calling the kdtree function
	double *query; //To store the random point for searching
	query = (double*) malloc(dim * sizeof(double));
	double *result_pt; //To store the nearest point to query
	result_pt = (double*) malloc(dim * sizeof(double));
	printf("KD Tree search results:\n");
	for(i=0;i<10;i++){ //Initializing the values to 0
		for(j=0;j<dim;j++){
			result_pt[j]=0.0;
			query[j] = rand() % 101; //Random generation of query point
		}
		printf("\tFor random point %d:\n",i+1);
		search_kdtree(dim,ndata,data,k,cluster_size,cluster_start,cluster_bdry,query,result_pt); //Call search function
	}
}

void kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid){
	int i=0,j=0;
	srand(time(NULL)); 
	for(i=0;i<dim*ndata;i++)
		data[i] = rand() % 101; //Randomly generate the points
	double min_val=0.0,max_val=0.0; //To store the min and max values for cluster boundaries
	//Initializing the values for initial cluster 0
	cluster_size[0] = ndata;
	cluster_start[0] = 0;
	for(i=0;i<dim;i++){
		min_val = max_val = data[cluster_start[0]+i];
		for(j=1;j<cluster_size[0];j++){
			if(data[i+j*dim]<min_val)
				min_val = data[i+j*dim];
			if(data[i+j*dim]>max_val)
				max_val = data[i+j*dim];
		}
		cluster_bdry[0][i] = min_val; 
		cluster_bdry[0][i+dim] = max_val;
	}
	for(cluster_counter=0;cluster_counter<k-1;cluster_counter++) // Calling the bipartition function for all the clusters
		bipartition(dim,cluster_start[cluster_counter],cluster_start[cluster_counter]+dim*cluster_size[cluster_counter]-1,data,cluster_size,cluster_start,cluster_bdry,cluster_centroid);
}

void bipartition(int dim, int i0, int im, double *data, int *cluster_size, int *cluster_start, double **cluster_bdry, double **cluster_centroid){
	int i=0,j=0;
	// Calculating the centroid for each cluster
	for(i=0;i<dim;i++){
		for(j=0;j<cluster_size[cluster_counter];j++)
			cluster_centroid[cluster_counter][i] = cluster_centroid[cluster_counter][i] + data[i0+i+j*dim];
		cluster_centroid[cluster_counter][i] = cluster_centroid[cluster_counter][i] / (double)cluster_size[cluster_counter];
	}
	double variance[dim]; //To store variance for each dimension
	for(i=0;i<dim;i++)
		variance[i]=0.0;
	// Variance calculation
	for(i=0;i<dim;i++){
		for(j=0;j<cluster_size[cluster_counter];j++)
			variance[i] = variance[i] + pow((data[i0+i+j*dim] - cluster_centroid[cluster_counter][i]),2);
		variance[i] = variance[i] / (double)cluster_size[cluster_counter];
	}
	int maxv_dim = 0; //To store the dimension having maximum variance
	double max_variance = variance[0]; //Initializing the maximum variance to first value in array
	// Finding the max variance
	for(i=1;i<dim;i++)
		if(variance[i]>max_variance)
			maxv_dim = i;
	double mean = cluster_centroid[cluster_counter][maxv_dim]; //Centroid value of the dimension having maximum variance
	int temp[cluster_size[cluster_counter]],temp1=0; //Temporary array to sort the data array using mean value
	for(i=0;i<cluster_size[cluster_counter];i++)
		temp[i] = -1;
	double temp2=0.0; //Used for swapping the value in data array
	for(i=0;i<cluster_size[cluster_counter];i++){
		if(data[i0+(maxv_dim)+i*dim]<mean)
			temp[i]=0; //Store 0 in temp array if data value is less than mean
		else
			temp[i]=1; //Store 1 in temp array if data value is greater than mean
	}
	//<-----------Sorting the data array using temp array
	int x=0,y=cluster_size[cluster_counter]-1;
	while(x<y){
		while(temp[x]==0)
			x++; //x keeps track of no. of zeroes in temp array
		while(temp[y]==1)
			y--;
		if(x<y){
			temp1=temp[x];
			temp[x]=temp[y];
			temp[y]=temp1;
			for(i=0;i<dim;i++){
				temp2=data[i0+i+x*dim];
				data[i0+i+x*dim]=data[i0+i+y*dim];
				data[i0+i+y*dim]=temp2;
			}
			x++;
			y--;
		}
	}
	// End of sorting the data array--------------->
	//Assigning the start and sizes of child clusters 
	cluster_size[2*cluster_counter+1] = x;
	cluster_size[2*cluster_counter+2] = cluster_size[cluster_counter] - x;
	cluster_start[2*cluster_counter+1] = i0;
	cluster_start[2*cluster_counter+2] = i0 + x*dim;
	//<---------- Boundary values calculation starts
	double min_val=0.0,max_val=0.0;
	for(i=0;i<dim;i++){
		// First child cluster boundaries calculation
		min_val = max_val = data[cluster_start[2*cluster_counter+1]+i];
		for(j=1;j<cluster_size[2*cluster_counter+1];j++){
			if(data[cluster_start[2*cluster_counter+1]+i+j*dim]<min_val)
				min_val = data[cluster_start[2*cluster_counter+1]+i+j*dim];
			if(data[cluster_start[2*cluster_counter+1]+i+j*dim]>max_val)
				max_val = data[cluster_start[2*cluster_counter+1]+i+j*dim];
		}
		cluster_bdry[2*cluster_counter+1][i] = min_val;
		cluster_bdry[2*cluster_counter+1][i+dim] = max_val;
		// Second child cluster boundaries calculation
		min_val = max_val = data[cluster_start[2*cluster_counter+2]+i];
		for(j=1;j<cluster_size[2*cluster_counter+2];j++){
			if(data[cluster_start[2*cluster_counter+2]+i+j*dim]<min_val)
				min_val = data[cluster_start[2*cluster_counter+2]+i+j*dim];
			if(data[cluster_start[2*cluster_counter+2]+i+j*dim]>max_val)
				max_val = data[cluster_start[2*cluster_counter+2]+i+j*dim];
		}
		cluster_bdry[2*cluster_counter+2][i] = min_val;
		cluster_bdry[2*cluster_counter+2][i+dim] = max_val;
	}
	// Boundary values calculation ends ------------->
}

void search_kdtree(int dim, int ndata, double *data, int k, int *cluster_size, int *cluster_start, double **cluster_bdry, double *query, double *result_pt){
	int i=0,j=0,visits=0,min_cluster=-1,cluster_visits=0;
	double cluster_min=-1.0,point_min=0.0,cluster_dist[2*k-1];
	for(i=0;i<2*k-1;i++)
		cluster_dist[i]=0.0; //Stores the distances between the query point and cluster boundaries
	// Cluster distances calculation 
	for(i=k-1;i<=2*k-2;i++){
		double cal_dist = 0.0;
		for(j=0;j<dim;j++){
			if(query[j]<cluster_bdry[i][j])
				cal_dist = cal_dist + pow((query[j]-cluster_bdry[i][j]),2);
			else if(query[j]>cluster_bdry[i][j+dim])
				cal_dist = cal_dist + pow((query[j]-cluster_bdry[i][j+dim]),2);
		}
		cluster_dist[i] = sqrt(cal_dist);
	}
	cluster_min=cluster_dist[k-1]; //Stores the min value of the cluster distances
	min_cluster = k-1; //Stores the index of the cluster having min distance from query point
	for(i=k;i<=2*k-2;i++){
		if(cluster_dist[i]<cluster_min){
			cluster_min = cluster_dist[i];
			min_cluster = i;
		}
	}
	//<---------- Calculating the distances from each point in the min cluster
	for(i=0;i<cluster_size[min_cluster];i++){
		double point_dist = 0.0;
		for(j=0;j<dim;j++)
			point_dist = point_dist + pow(query[j]-data[cluster_start[min_cluster]+i*dim+j],2);
		point_dist = sqrt(point_dist);
		if(i==0)
			point_min = point_dist;
		if(point_dist<point_min){
			point_min = point_dist; //point_min has the minimum distance between the query point and one of the point in min cluster
			for(j=0;j<dim;j++)
				result_pt[j] = data[cluster_start[min_cluster]+i*dim+j]; //Stores the point having the min distance from query point
		}
		visits++; //Keeps track of the total number of points visited
	}
	cluster_visits++; //Keeps track of the total number of clusters visited
	//<---------- Finding whether the min distance exists in any other cluster
	for(i=k-1;i<=2*k-2;i++){
		if(i==min_cluster)
			continue;
		if(cluster_dist[i]<=point_min){
			for(j=0;j<cluster_size[i];j++){
				double new_dist = 0.0;
				int z=0;
				for(z=0;z<dim;z++)
					new_dist = new_dist + pow(query[z]-data[cluster_start[i]+j*dim+z],2);
				if(sqrt(new_dist)<point_min){
					point_min = sqrt(new_dist);
					for(z=0;z<dim;z++)
						result_pt[z] = data[cluster_start[i]+j*dim+z];
				}
				visits++;
			}
			cluster_visits++;
		}
	}
	// End of finding the min distance. Final output is obtained here --------->
	printf("\t\tTotal Points Visited: %d\n\t\tMinimum Distance: %f\n",visits,point_min);
}