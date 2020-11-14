function keep_guis
fig_h = permute( findobj( 0, 'Type', 'Figure' ), [2,1] );
for fh = fig_h
    uih = findobj( fh, 'Type', 'uicontrol' );
    if isempty( uih )
        delete( fh );
    end
end
end