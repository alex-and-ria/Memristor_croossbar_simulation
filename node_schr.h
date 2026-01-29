









































void node_analyzer(unsigned int* row, unsigned int *col, unsigned int *len, unsigned int* nds_td, unsigned int* nds_n, double* nb_koef, unsigned int* n_indp);
void star_mesh_base(unsigned int *row,unsigned int** rw, unsigned int *col,unsigned int** cl, double *val,double** vl, unsigned int *len,unsigned int *ln, unsigned int *nds_td, unsigned int *nds_n);
void dense_rdct(unsigned int *row,unsigned int*** rw_, unsigned int *col,unsigned int*** cl_, double *val,double*** vl_, unsigned int *len,unsigned int **ln_, unsigned int *nds_td, unsigned int *nds_n,double th_nb_koef, unsigned int *nsd_td1, unsigned int nds_n1);
int tst_pnt(double** adj_m,unsigned int* dim);
