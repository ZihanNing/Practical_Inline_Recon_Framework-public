function  res = ctranspose(obj)

obj.adjoint = xor(obj.adjoint,1);
res = obj;
