function rsnr = calcRSNR(x,x0)
    rsnr = 20*log10(norm(x0(:))/norm(x0(:)-x(:)));
end