.arch armv8-a
.global fmla_test
.type fmla_test %function
fmla_test:
	//Param# Register Datatype Name
	// 1     X0       long     k
	k_iter .req x0
	mov x1, #1
	mov x2, #1

	fmov d0, #1.0 // set all vector registers to 1.0
	dup v0.2d, v0.2d[0]
	dup v1.2d, v0.2d[0]
	dup v2.2d, v0.2d[0]
	dup v3.2d, v0.2d[0]
	dup v4.2d, v0.2d[0]
	dup v5.2d, v0.2d[0]
	dup v6.2d, v0.2d[0]
	dup v7.2d, v0.2d[0]
	dup v8.2d, v0.2d[0]
	dup v9.2d, v0.2d[0]
	dup v10.2d, v0.2d[0]
	dup v11.2d, v0.2d[0]
	dup v12.2d, v0.2d[0]
	dup v13.2d, v0.2d[0]
	dup v14.2d, v0.2d[0]
	dup v15.2d, v0.2d[0]
	dup v16.2d, v0.2d[0]
	dup v17.2d, v0.2d[0]
	dup v18.2d, v0.2d[0]
	dup v19.2d, v0.2d[0]
	dup v20.2d, v0.2d[0]
	dup v21.2d, v0.2d[0]
	dup v22.2d, v0.2d[0]
	dup v23.2d, v0.2d[0]
	dup v24.2d, v0.2d[0]
	dup v25.2d, v0.2d[0]
	dup v26.2d, v0.2d[0]
	dup v27.2d, v0.2d[0]
	dup v28.2d, v0.2d[0]
	dup v29.2d, v0.2d[0]
	dup v30.2d, v0.2d[0]
	dup v31.2d, v0.2d[0]
testit:
	sub k_iter, k_iter, #1
	fmla v2.2d, v0.2d, v1.2d
	add x3, x26, x27
	fmla v3.2d, v0.2d, v1.2d
	add x4, x26, x27
	fmla v4.2d, v0.2d, v1.2d
	add x5, x26, x27
	fmla v5.2d, v0.2d, v1.2d
	add x6, x26, x27
	fmla v6.2d, v0.2d, v1.2d
	add x7, x26, x27
	fmla v7.2d, v0.2d, v1.2d
	add x8, x26, x27
	fmla v8.2d, v0.2d, v1.2d
	add x9, x26, x27
	fmla v9.2d, v0.2d, v1.2d
	add x10, x26, x27

	fmla v10.2d, v0.2d, v1.2d
	add x11, x26, x27
	fmla v11.2d, v0.2d, v1.2d
	add x12, x26, x27
	fmla v12.2d, v0.2d, v1.2d
	add x13, x26, x27
	fmla v13.2d, v0.2d, v1.2d
	add x14, x26, x27
	fmla v14.2d, v0.2d, v1.2d
	add x15, x26, x27
	fmla v15.2d, v0.2d, v1.2d
	add x16, x26, x27
	fmla v16.2d, v0.2d, v1.2d
	add x17, x26, x27
	fmla v17.2d, v0.2d, v1.2d
	add x18, x26, x27

	cbnz k_iter, testit


	ret

