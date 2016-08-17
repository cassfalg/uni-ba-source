.arch armv8-a
.global my_dgemm_simple_4x4_cm
.type my_dgemm_simple_4x4_cm %function
my_dgemm_simple_4x4_cm:
	//Param# Register Datatype Name
	// 1     X0       long     k
	// 2     D0       double   alpha
	// 3     X1       double   *A
	// 4     X2       double   *B
	// 5     D1       double   beta
	// 6     X3       double   *C
	// 7     X4       long     incRowC
	// 8     X5       long     incColC

	/*
	Preliminaries:
	save alpha and beta
	define register aliases for a, b and C_
	zero C_
	*/
	ins v0.d[1], v0.d[0] // copy alpha to the higher lane of v0
	ins v0.d[0], v1.d[0] // copy beta to the lower lane of v0
	alpha_beta .req v0 // alpha in higher lane, beta in lower lane
	beta .req d0 // same as v0.d[0] for reading
	a1 .req v1
	qa1 .req q1 // need q name alias for load/store, v for data processing
	a2 .req v2
	qa2 .req q2
	b_1 .req v3 // b1 and b2 are built in registers, can't be redefined
	qb_1 .req q3
	b_2 .req v4
	qb_2 .req q4
	c_1_1 .req v24
	c_2_1 .req v25
	c_3_1 .req v26
	c_4_1 .req v27
	c_1_2 .req v28
	c_2_2 .req v29
	c_3_2 .req v30
	c_4_2 .req v31
	dup c_1_1.2d, xzr // initialize C_ with zeros
	dup c_2_1.2d, xzr
	dup c_3_1.2d, xzr
	dup c_4_1.2d, xzr
	dup c_1_2.2d, xzr
	dup c_2_2.2d, xzr
	dup c_3_2.2d, xzr
	dup c_4_2.2d, xzr

	/*
	Accumulate - loop over k:
	load current row vector a_k into a
	load current column vector b_k into b
	calculate the dyadic product
	*/
	k_iter .req x6
	add k_iter, x0, #1
	Ap .req x7
	mov Ap, x1
	Bp .req x8
	mov Bp, x2
accumulate:
	sub k_iter, k_iter, #1
	ldr qa1, [Ap], #16 // load a and b
	ldr qa2, [Ap], #16
	ldr qb_1, [Bp], #16
	ldr qb_2, [Bp], #16
	fmla c_1_1.2d, b_1.2d, a1.2d[0] // calculate dyadic product
	fmla c_2_1.2d, b_1.2d, a1.2d[1]
	fmla c_3_1.2d, b_1.2d, a2.2d[0]
	fmla c_4_1.2d, b_1.2d, a2.2d[1]
	fmla c_1_2.2d, b_2.2d, a1.2d[0]
	fmla c_2_2.2d, b_2.2d, a1.2d[1]
	fmla c_3_2.2d, b_2.2d, a2.2d[0]
	fmla c_4_2.2d, b_2.2d, a2.2d[1]
	cbnz k_iter, accumulate

	/*
	Scale results and write back:
	load C from memory
	scale C by beta
	scale C_ by alpha and add to C
	store C to memory
	*/
	Cp .req x3
	c1_1 .req v16 // registers used to store C from memory
	c2_1 .req v17
	c3_1 .req v18
	c4_1 .req v19
	c1_2 .req v20
	c2_2 .req v21
	c3_2 .req v22
	c4_2 .req v23
	incRowCb .req x13
	lsl incRowCb, x4, #3 // incRowC * sizeof(double) (we need it in bytes)
	incColCb .req x14
	lsl incColCb, x5, #3 // incColC * sizeof(double)
	Cp1 .req x9 // pointer on column of C
	add Cp1, Cp, xzr // adress column1
	Cp2 .req x10
	add Cp2, Cp1, incColCb // adress column2
	Cp3 .req x11
	add Cp3, Cp2, incColCb // adress column3
	Cp4 .req x12
	add Cp4, Cp3, incColCb // adress column4

	/*
	check if C needs to be initialized with zero or loaded from memory and scaled
	*/
	fcmp beta, #0.0
	b.eq beta_is_zero

	ld1 {c1_1.d}[0], [cp1], incRowCb // load column 1
	ld1 {c2_1.d}[0], [cp1], incRowCb
	ld1 {c3_1.d}[0], [cp1], incRowCb
	ld1 {c4_1.d}[0], [cp1], incRowCb
	ld1 {c1_1.d}[1], [cp2], incRowCb // load column 2
	ld1 {c2_1.d}[1], [cp2], incRowCb
	ld1 {c3_1.d}[1], [cp2], incRowCb
	ld1 {c4_1.d}[1], [cp2], incRowCb
	ld1 {c1_2.d}[0], [cp3], incRowCb // load column 3
	ld1 {c2_2.d}[0], [cp3], incRowCb
	ld1 {c3_2.d}[0], [cp3], incRowCb
	ld1 {c4_2.d}[0], [cp3], incRowCb
	ld1 {c1_2.d}[1], [cp4], incRowCb // load column 4
	ld1 {c2_2.d}[1], [cp4], incRowCb
	ld1 {c3_2.d}[1], [cp4], incRowCb
	ld1 {c4_2.d}[1], [cp4], incRowCb
	fmul c1_1.2d, c1_1.2d, alpha_beta.2d[0] // scale C by beta
	fmul c2_1.2d, c2_1.2d, alpha_beta.2d[0]
	fmul c3_1.2d, c3_1.2d, alpha_beta.2d[0]
	fmul c4_1.2d, c4_1.2d, alpha_beta.2d[0]
	fmul c1_2.2d, c1_2.2d, alpha_beta.2d[0]
	fmul c2_2.2d, c2_2.2d, alpha_beta.2d[0]
	fmul c3_2.2d, c3_2.2d, alpha_beta.2d[0]
	fmul c4_2.2d, c4_2.2d, alpha_beta.2d[0]
	b add_alpha_times_C_ // skip beta=0 case and continue, adding alpha * C_
beta_is_zero:
	dup c1_1.2d, xzr // initialize C with zeros
	dup c2_1.2d, xzr
	dup c3_1.2d, xzr
	dup c4_1.2d, xzr
	dup c1_2.2d, xzr
	dup c2_2.2d, xzr
	dup c3_2.2d, xzr
	dup c4_2.2d, xzr
add_alpha_times_C_:
	fmla c1_1.2d, c_1_1.2d, alpha_beta.2d[1] // C = C + alpha * C_
	fmla c2_1.2d, c_2_1.2d, alpha_beta.2d[1]
	fmla c3_1.2d, c_3_1.2d, alpha_beta.2d[1]
	fmla c4_1.2d, c_4_1.2d, alpha_beta.2d[1]
	fmla c1_2.2d, c_1_2.2d, alpha_beta.2d[1]
	fmla c2_2.2d, c_2_2.2d, alpha_beta.2d[1]
	fmla c3_2.2d, c_3_2.2d, alpha_beta.2d[1]
	fmla c4_2.2d, c_4_2.2d, alpha_beta.2d[1]
	add Cp1, Cp, xzr // adress column1
	add Cp2, Cp1, incColCb // adress column2
	add Cp3, Cp2, incColCb // adress column3
	add Cp4, Cp3, incColCb // adress column4
	st1 {c1_1.d}[0], [cp1], incRowCb // store column 1
	st1 {c2_1.d}[0], [cp1], incRowCb
	st1 {c3_1.d}[0], [cp1], incRowCb
	st1 {c4_1.d}[0], [cp1], incRowCb
	st1 {c1_1.d}[1], [cp2], incRowCb // store column 2
	st1 {c2_1.d}[1], [cp2], incRowCb
	st1 {c3_1.d}[1], [cp2], incRowCb
	st1 {c4_1.d}[1], [cp2], incRowCb
	st1 {c1_2.d}[0], [cp3], incRowCb // store column 3
	st1 {c2_2.d}[0], [cp3], incRowCb
	st1 {c3_2.d}[0], [cp3], incRowCb
	st1 {c4_2.d}[0], [cp3], incRowCb
	st1 {c1_2.d}[1], [cp4], incRowCb // store column 4
	st1 {c2_2.d}[1], [cp4], incRowCb
	st1 {c3_2.d}[1], [cp4], incRowCb
	st1 {c4_2.d}[1], [cp4], incRowCb

	ret

