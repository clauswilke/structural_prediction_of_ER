	function select(k,n,arr)
	 integer k,n
	 double precision select,arr(n)
!	 Returns the kth smallest value in the array arr(1:n). The input array will be rearranged to have this value in location arr(k), 
!	 with all smaller elements move to arr(1:k-1) (in arbitrary order) and all larger elements in arr(K+1:n) (also in arbitrary order)
	 integer i, ir, j, l, mid
	 double precision a, temp
	 l = 1
	 ir = n
1	 if (ir-l .le. 1) then
		if (ir-l .eq. 1) then
			if (arr(ir) .lt. arr(l)) then
				temp = arr(l)
				arr(l) = arr(ir)
				arr(ir) = temp
			end if
		end if
		select = arr(k)
		return
	 else
		mid = (l + ir)/2
		temp = arr(mid)
		arr(mid) = arr(l+1)
		arr(l+1) = temp
		if (arr(l+1) .gt. arr(ir)) then
			temp = arr(l+1)
			arr(l+1) = arr(ir)
			arr(ir) = temp
		end if
		if (arr(l) .gt. arr(ir)) then
			temp = arr(l)
			arr(l) = arr(ir)
			arr(ir) = temp
		end if
		if (arr(l+1) .gt. arr(l)) then
			temp = arr(l+1)
			arr(l+1) = arr(l)
			arr(l) = temp
		end if
		i = l + 1
		j = ir
		a = arr(l)
3		continue
			i = i + 1
		if (arr(i) .lt. a) goto 3
4		continue
			j = j - 1
		if (arr(j) .gt. a) goto 4
		if (j .lt. i) goto 5
		temp = arr(i)
		arr(i) = arr(j)
		arr(j) = temp
		goto 3
5		arr(l) = arr(j)
		arr(j) = a
		if (j .ge. k) ir=j-1
		if (j .le. k) l=i
	 end if
	 goto 1
	end
		
	
	 
