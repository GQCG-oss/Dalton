
/* routines for generation of Lebedev grid */
void ld0026_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0038_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0050_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0074_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0086_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0110_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0146_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0170_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0194_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0230_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0266_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0302_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0350_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0590_(real* x, real* y, real* z, real* w, integer *pnt);
void ld0974_(real* x, real* y, real* z, real* w, integer *pnt);
void ld1202_(real* x, real* y, real* z, real* w, integer *pnt);
void ld1454_(real* x, real* y, real* z, real* w, integer *pnt);

struct leb_gen_{
  integer point_cnt;
  void (*func)(real* x, real* y, real* z, real* w, integer *pnt);
};

struct leb_gen_ leb_gen[] = {
 {  26, ld0026_},
 {  38, ld0038_},
 {  50, ld0050_},
 {  74, ld0074_},
 {  86, ld0086_},
 { 110, ld0110_},
 { 146, ld0146_},
 { 170, ld0170_},
 { 194, ld0194_},
 { 230, ld0230_},
 { 266, ld0266_},
 { 302, ld0302_},
 { 350, ld0350_},
 { 590, ld0590_},
 { 974, ld0974_},
 {1202, ld1202_},
 {1454, ld1454_},
};  
