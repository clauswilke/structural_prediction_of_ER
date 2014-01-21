! number of space-separated items in a line
integer function nitems(line)
   character line*(*)    
   logical back
   integer length
   
   back = .true.        
   length = len_trim(line)    
   k = index(line(1:length), ' ', back)
   if (k == 0) then
       nitems = 0
       return
   end if    
   
   nitems = 1
   do 
       ! starting with the right most blank space, 
       ! look for the next non-space character down
       ! indicating there is another item in the line
       do
           if (k <= 0) exit
           
           if (line(k:k) == ' ') then
               k = k - 1
               cycle
           else
               nitems = nitems + 1
               exit
           end if
           
       end do
       
       ! once a non-space character is found,
       ! skip all adjacent non-space character
       do
           if ( k<=0 ) exit
           
           if (line(k:k) /= ' ') then
               k = k - 1
               cycle
           end if
           
           exit
           
       end do
       
       if (k <= 0) exit
           
   end do
end function nitems 
