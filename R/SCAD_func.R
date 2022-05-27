SCAD_func<-function(yvpen,lambda_ch=0.01,gamma=3)
{
  if(gamma<=2){gamma = 2.01;print("Gamma needs > 2!!!")}
  ynew = sign(yvpen)*(abs(yvpen)-lambda_ch)*(abs(yvpen)>=lambda_ch)*(abs(yvpen)<=2*lambda_ch)+ 
    yvpen*(abs(yvpen) > gamma*lambda_ch) + 
    ((gamma-1)*yvpen-sign(yvpen)*gamma*lambda_ch)/(gamma-2)*(abs(yvpen)<=gamma*lambda_ch)*(abs(yvpen)>2*lambda_ch)
  if(sd(ynew) < 0.0000001){print("Parmeters are not correctly specified!");return(ynew)}
  return(ynew)#/sd(ynew)*sd(yvpen))
}
