function h = display_shape(shape, colorMap)

  if nargin < 2
    h = trisurf(shape.surface.TRIV, shape.surface.VERT(:,1), shape.surface.VERT(:,2), shape.surface.VERT(:,3), 'FaceColor', 'interp');
  else
    h = trisurf(shape.surface.TRIV, shape.surface.VERT(:,1), shape.surface.VERT(:,2), shape.surface.VERT(:,3), colorMap,'FaceColor', 'interp');

  end
  
  set(h, 'edgecolor', 'none');
  set(h, 'FaceAlpha', 0.5);
  axis equal; axis off; hold on;

end