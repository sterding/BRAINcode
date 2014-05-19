# for median
function median(v) 
{ 
    c=asort(v,j); 
    if (c % 2) return j[(c+1)/2]; 
    else return (j[c/2+1]+j[c/2])/2; 
}

# for trimmed mean  
function trimmedMean(v, p) 
{ 
    c=asort(v,j); 
    a=int(c*p);
    s=0;
    for(i=a+1;i<=(c-a);i++) s+=j[i];
    return s/(c-2*a); 
} 
