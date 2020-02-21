% computes the error using the ki
function [masserr, engerr] = compute_energy_errors(q,qe,q_ref,grid_obj)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    masse=0.0; massn=0.0;
    enge=0.0; engn=0.0;
    for ie=1:nel
        q1 = grid_obj.normvdm*q(1,:,ie)';
        q2 = grid_obj.normvdm*q(2,:,ie)';
        qr = grid_obj.normvdm*q_ref(1,:,ie)';
        qe1 = grid_obj.normvdm*qe(1,:,ie)';
        qe2 = grid_obj.normvdm*qe(2,:,ie)';
        for i=1:2*ngl+1
            x = grid_obj.normgrd(i,ie);
            wq = grid_obj.normwq(i,ie);

            masse = masse + wq*( qe1(i) );
            massn = massn + wq*( q1(i) );
            
            enge = enge + wq*0.5*( qe1(i)^2 + safe_div(qe2(i)^2, qe1(i), 0.0) + qe1(i)*qr(i));
            engn = engn + wq*0.5*( q1(i)^2 + safe_div(q2(i)^2, q1(i), 0.0) + q1(i)*qr(i) );
            
        end
    end
    masserr = (massn - masse)/masse;
    engerr = (engn - enge)/enge;
end