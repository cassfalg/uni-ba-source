.arch armv8-a
.global my_dgemm_4x8ugemm
.type my_dgemm_4x8ugemm %function
my_dgemm_4x8ugemm:
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
	push registers v8 to v15 (callee saved)
	save alpha and beta
	define register aliases for a, b and C_
	zero C_
	*/
	stp d8, d9, [sp, #-16]! // push registers to stack
	stp d10, d11, [sp, #-16]!
	stp d12, d13, [sp, #-16]!
	stp d14, d15, [sp, #-16]!
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
	b_3 .req v5
	qb_3 .req q5
	b_4 .req v6
	qb_4 .req q6
	c_1_1 .req v16
	c_2_1 .req v17
	c_1_2 .req v18
	c_2_2 .req v19
	c_1_3 .req v20
	c_2_3 .req v21
	c_1_4 .req v22
	c_2_4 .req v23
	c_1_5 .req v24
	c_2_5 .req v25
	c_1_6 .req v26
	c_2_6 .req v27
	c_1_7 .req v28
	c_2_7 .req v29
	c_1_8 .req v30
	c_2_8 .req v31
	dup c_1_1.2d, xzr // initialize C_ with zeros
	dup c_2_1.2d, xzr
	dup c_1_2.2d, xzr
	dup c_2_2.2d, xzr
	dup c_1_3.2d, xzr
	dup c_2_3.2d, xzr
	dup c_1_4.2d, xzr
	dup c_2_4.2d, xzr
	dup c_1_5.2d, xzr
	dup c_2_5.2d, xzr
	dup c_1_6.2d, xzr
	dup c_2_6.2d, xzr
	dup c_1_7.2d, xzr
	dup c_2_7.2d, xzr
	dup c_1_8.2d, xzr
	dup c_2_8.2d, xzr

	/*
	Accumulate - loop over k:
	load current row vector a_k into a
	load current column vector b_k into b
	calculate the dyadic product
	*/
	k_iter .req x6
	add k_iter, x0, #0
	Ap .req x7
	mov Ap, x1
	Bp .req x9
	mov Bp, x2
accumulate:
	sub k_iter, k_iter, #1
	ldr qa1, [Ap], #16 // load a and b
	ldr qa2, [Ap], #16
	ldr qb_1, [Bp], #16
	ldr qb_2, [Bp], #16
	ldr qb_3, [Bp], #16
	ldr qb_4, [Bp], #16
	fmla c_1_1.2d, a1.2d, b_1.2d[0] // calculate dyadic product
	fmla c_2_1.2d, a2.2d, b_1.2d[0]
	fmla c_1_2.2d, a1.2d, b_1.2d[1]
	fmla c_2_2.2d, a2.2d, b_1.2d[1]
	fmla c_1_3.2d, a1.2d, b_2.2d[0]
	fmla c_2_3.2d, a2.2d, b_2.2d[0]
	fmla c_1_4.2d, a1.2d, b_2.2d[1]
	fmla c_2_4.2d, a2.2d, b_2.2d[1]
	fmla c_1_5.2d, a1.2d, b_3.2d[0]
	fmla c_2_5.2d, a2.2d, b_3.2d[0]
	fmla c_1_6.2d, a1.2d, b_3.2d[1]
	fmla c_2_6.2d, a2.2d, b_3.2d[1]
	fmla c_1_7.2d, a1.2d, b_4.2d[0]
	fmla c_2_7.2d, a2.2d, b_4.2d[0]
	fmla c_1_8.2d, a1.2d, b_4.2d[1]
	fmla c_2_8.2d, a2.2d, b_4.2d[1]
	cbnz k_iter, accumulate

	/*
	Scale results and write back:
	load C from memory
	scale C by beta
	scale C_ by alpha and add to C
	store C to memory
	*/
	Cp .req x3
	c1_1 .req v7 // registers used to store C from memory
	c2_1 .req v8
	c1_2 .req v9
	c2_2 .req v10
	c1_3 .req v11
	c2_3 .req v12
	c1_4 .req v13
	c2_4 .req v14
	incRowCb .req x10
	lsl incRowCb, x4, #3 // incRowC * sizeof(double) (we need it in bytes)
	incColCb .req x11
	lsl incColCb, x5, #3 // incColC * sizeof(double)
	Cp1 .req x12 // pointer on column of C
	add Cp1, Cp, xzr // adress column1
	Cp2 .req x13
	add Cp2, Cp1, incColCb // adress column2
	Cp3 .req x14
	add Cp3, Cp2, incColCb // adress column3
	Cp4 .req x15 //use x6 because it is free and we need Cp4 twice for the second half of C/C_
	add Cp4, Cp3, incColCb // adress column4
	Cp4_sav .req x7
	mov Cp4_sav, Cp4 // we need this for the second half of C/C_

	/*
	check if C needs to be initialized with zero or loaded from memory and scaled
	*/
	fcmp beta, #0.0
	b.eq beta_is_zero_1

	ld1 {c1_1.d}[0], [cp1], incRowCb // load column 1
	ld1 {c1_1.d}[1], [cp1], incRowCb
	ld1 {c2_1.d}[0], [cp1], incRowCb
	ld1 {c2_1.d}[1], [cp1], incRowCb
	ld1 {c1_2.d}[0], [cp2], incRowCb // load column 2
	ld1 {c1_2.d}[1], [cp2], incRowCb
	ld1 {c2_2.d}[0], [cp2], incRowCb
	ld1 {c2_2.d}[1], [cp2], incRowCb
	ld1 {c1_3.d}[0], [cp3], incRowCb // load column 3
	ld1 {c1_3.d}[1], [cp3], incRowCb
	ld1 {c2_3.d}[0], [cp3], incRowCb
	ld1 {c2_3.d}[1], [cp3], incRowCb
	ld1 {c1_4.d}[0], [cp4], incRowCb // load column 4
	ld1 {c1_4.d}[1], [cp4], incRowCb
	ld1 {c2_4.d}[0], [cp4], incRowCb
	ld1 {c2_4.d}[1], [cp4], incRowCb
	fmul c1_1.2d, c1_1.2d, alpha_beta.2d[0] // scale C by beta
	fmul c2_1.2d, c2_1.2d, alpha_beta.2d[0]
	fmul c1_2.2d, c1_2.2d, alpha_beta.2d[0]
	fmul c2_2.2d, c2_2.2d, alpha_beta.2d[0]
	fmul c1_3.2d, c1_3.2d, alpha_beta.2d[0]
	fmul c2_3.2d, c2_3.2d, alpha_beta.2d[0]
	fmul c1_4.2d, c1_4.2d, alpha_beta.2d[0]
	fmul c2_4.2d, c2_4.2d, alpha_beta.2d[0]
	b add_alpha_times_C_1 // skip beta=0 case and continue, adding alpha * C_
beta_is_zero_1:
	dup c1_1.2d, xzr // initialize C with zeros
	dup c2_1.2d, xzr
	dup c1_2.2d, xzr
	dup c2_2.2d, xzr
	dup c1_3.2d, xzr
	dup c2_3.2d, xzr
	dup c1_4.2d, xzr
	dup c2_4.2d, xzr
add_alpha_times_C_1:
	fmla c1_1.2d, c_1_1.2d, alpha_beta.2d[1] // C = C + alpha * C_
	fmla c2_1.2d, c_2_1.2d, alpha_beta.2d[1]
	fmla c1_2.2d, c_1_2.2d, alpha_beta.2d[1]
	fmla c2_2.2d, c_2_2.2d, alpha_beta.2d[1]
	fmla c1_3.2d, c_1_3.2d, alpha_beta.2d[1]
	fmla c2_3.2d, c_2_3.2d, alpha_beta.2d[1]
	fmla c1_4.2d, c_1_4.2d, alpha_beta.2d[1]
	fmla c2_4.2d, c_2_4.2d, alpha_beta.2d[1]
	add Cp1, Cp, xzr // adress column1
	add Cp2, Cp1, incColCb // adress column2
	add Cp3, Cp2, incColCb // adress column3
	add Cp4, Cp3, incColCb // adress column4
	st1 {c1_1.d}[0], [cp1], incRowCb // store column 1
	st1 {c1_1.d}[1], [cp1], incRowCb
	st1 {c2_1.d}[0], [cp1], incRowCb
	st1 {c2_1.d}[1], [cp1], incRowCb
	st1 {c1_2.d}[0], [cp2], incRowCb // store column 2
	st1 {c1_2.d}[1], [cp2], incRowCb
	st1 {c2_2.d}[0], [cp2], incRowCb
	st1 {c2_2.d}[1], [cp2], incRowCb
	st1 {c1_3.d}[0], [cp3], incRowCb // store column 3
	st1 {c1_3.d}[1], [cp3], incRowCb
	st1 {c2_3.d}[0], [cp3], incRowCb
	st1 {c2_3.d}[1], [cp3], incRowCb
	st1 {c1_4.d}[0], [cp4], incRowCb // store column 4
	st1 {c1_4.d}[1], [cp4], incRowCb
	st1 {c2_4.d}[0], [cp4], incRowCb
	st1 {c2_4.d}[1], [cp4], incRowCb


	/*
	Now do the same for the second part of C_
	*/
	//Cp .req x3
	c1_5 .req v7 // registers used to store C from memory
	c2_5 .req v8
	c1_6 .req v9
	c2_6 .req v10
	c1_7 .req v11
	c2_7 .req v12
	c1_8 .req v13
	c2_8 .req v14
	//incRowCb .req x10
	//lsl incRowCb, x4, #3 // incRowC * sizeof(double) (we need it in bytes)
	//incColCb .req x11
	//lsl incColCb, x5, #3 // incColC * sizeof(double)
	Cp5 .req x12
	add Cp5, Cp4_sav, incColCb // adress column5, reuse Cp4 from before. Sames as Cp8 later
	Cp6 .req x13
	add Cp6, Cp5, incColCb // adress column6
	Cp7 .req x14
	add Cp7, Cp6, incColCb // adress column7
	Cp8 .req x15
	add Cp8, Cp7, incColCb // adress column8

	/*
	check if C needs to be initialized with zero or loaded from memory and scaled
	*/
	fcmp beta, #0.0
	b.eq beta_is_zero_2

	ld1 {c1_5.d}[0], [cp5], incRowCb // load column 1
	ld1 {c1_5.d}[1], [cp5], incRowCb
	ld1 {c2_5.d}[0], [cp5], incRowCb
	ld1 {c2_5.d}[1], [cp5], incRowCb
	ld1 {c1_6.d}[0], [cp6], incRowCb // load column 2
	ld1 {c1_6.d}[1], [cp6], incRowCb
	ld1 {c2_6.d}[0], [cp6], incRowCb
	ld1 {c2_6.d}[1], [cp6], incRowCb
	ld1 {c1_7.d}[0], [cp7], incRowCb // load column 3
	ld1 {c1_7.d}[1], [cp7], incRowCb
	ld1 {c2_7.d}[0], [cp7], incRowCb
	ld1 {c2_7.d}[1], [cp7], incRowCb
	ld1 {c1_8.d}[0], [cp8], incRowCb // load column 4
	ld1 {c1_8.d}[1], [cp8], incRowCb
	ld1 {c2_8.d}[0], [cp8], incRowCb
	ld1 {c2_8.d}[1], [cp8], incRowCb
	fmul c1_5.2d, c1_5.2d, alpha_beta.2d[0] // scale C by beta
	fmul c2_5.2d, c2_5.2d, alpha_beta.2d[0]
	fmul c1_6.2d, c1_6.2d, alpha_beta.2d[0]
	fmul c2_6.2d, c2_6.2d, alpha_beta.2d[0]
	fmul c1_7.2d, c1_7.2d, alpha_beta.2d[0]
	fmul c2_7.2d, c2_7.2d, alpha_beta.2d[0]
	fmul c1_8.2d, c1_8.2d, alpha_beta.2d[0]
	fmul c2_8.2d, c2_8.2d, alpha_beta.2d[0]
	b add_alpha_times_C_2 // skip beta=0 case and continue, adding alpha * C_
beta_is_zero_2:
	dup c1_5.2d, xzr // initialize C with zeros
	dup c2_5.2d, xzr
	dup c1_6.2d, xzr
	dup c2_6.2d, xzr
	dup c1_7.2d, xzr
	dup c2_7.2d, xzr
	dup c1_8.2d, xzr
	dup c2_8.2d, xzr
add_alpha_times_C_2:
	fmla c1_5.2d, c_1_5.2d, alpha_beta.2d[1] // C = C + alpha * C_
	fmla c2_5.2d, c_2_5.2d, alpha_beta.2d[1]
	fmla c1_6.2d, c_1_6.2d, alpha_beta.2d[1]
	fmla c2_6.2d, c_2_6.2d, alpha_beta.2d[1]
	fmla c1_7.2d, c_1_7.2d, alpha_beta.2d[1]
	fmla c2_7.2d, c_2_7.2d, alpha_beta.2d[1]
	fmla c1_8.2d, c_1_8.2d, alpha_beta.2d[1]
	fmla c2_8.2d, c_2_8.2d, alpha_beta.2d[1]
	add Cp5, Cp4_sav, incColCb // adress column5
	add Cp6, Cp5, incColCb // adress column6
	add Cp7, Cp6, incColCb // adress column7
	add Cp8, Cp7, incColCb // adress column8
	st1 {c1_5.d}[0], [cp5], incRowCb // store column 5
	st1 {c1_5.d}[1], [cp5], incRowCb
	st1 {c2_5.d}[0], [cp5], incRowCb
	st1 {c2_5.d}[1], [cp5], incRowCb
	st1 {c1_6.d}[0], [cp6], incRowCb // store column 6
	st1 {c1_6.d}[1], [cp6], incRowCb
	st1 {c2_6.d}[0], [cp6], incRowCb
	st1 {c2_6.d}[1], [cp6], incRowCb
	st1 {c1_7.d}[0], [cp7], incRowCb // store column 7
	st1 {c1_7.d}[1], [cp7], incRowCb
	st1 {c2_7.d}[0], [cp7], incRowCb
	st1 {c2_7.d}[1], [cp7], incRowCb
	st1 {c1_8.d}[0], [cp8], incRowCb // store column 8
	st1 {c1_8.d}[1], [cp8], incRowCb
	st1 {c2_8.d}[0], [cp8], incRowCb
	st1 {c2_8.d}[1], [cp8], incRowCb

	ldp d14, d15, [sp], #16 // pop registers from stack
	ldp d12, d13, [sp], #16
	ldp d10, d11, [sp], #16
	ldp d8, d9, [sp], #16

	ret

