









































void node_analyzer(unsigned int* row, unsigned int *col, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef, unsigned int* n_indp);
void star_mesh_base(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n);
void dense_rdct(unsigned int *row, unsigned long long int* rw_v, unsigned int *col, unsigned long long int* cl_v, double *val, unsigned long long int* vl_v, unsigned int *len,unsigned int **ln_, unsigned int *nds_td, unsigned int *nds_n,double th_nb_koef, unsigned int *nsd_td1, unsigned int nds_n1,int mode_dbg, unsigned int num_iter);
void get_dbg_arr(unsigned long long int *rw_v,unsigned long long int* cl_v, unsigned long long int* vl_v, unsigned int** rw0, unsigned int** cl0, double **vl0);
//void get_pointer(long long unsigned int* ptr_tmp, unsigned int num_out);
//int tst_pnt(double** adj_m,unsigned int* dim);
int tst_pnt(unsigned int dim);
