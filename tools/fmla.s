.arch armv8-a
.global fmla_test
.type fmla_test %function
fmla_test:
	//Param# Register Datatype Name
	// 1     X0       long     k
	k_iter .req x0

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
	fmla v3.2d, v0.2d, v1.2d
	fmla v4.2d, v0.2d, v1.2d
	fmla v5.2d, v0.2d, v1.2d
	fmla v6.2d, v0.2d, v1.2d
	fmla v7.2d, v0.2d, v1.2d
	fmla v8.2d, v0.2d, v1.2d
	fmla v9.2d, v0.2d, v1.2d

	fmla v10.2d, v0.2d, v1.2d
	fmla v11.2d, v0.2d, v1.2d
	fmla v12.2d, v0.2d, v1.2d
	fmla v13.2d, v0.2d, v1.2d
	fmla v14.2d, v0.2d, v1.2d
	fmla v15.2d, v0.2d, v1.2d
	fmla v16.2d, v0.2d, v1.2d
	fmla v17.2d, v0.2d, v1.2d

	cbnz k_iter, testit


	ret

