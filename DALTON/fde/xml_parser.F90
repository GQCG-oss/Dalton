!  Copyright (C) 2018 Andre Severo Pereira Gomes, Christoph Jacob, Lucas Visscher and collaborators
!
!  This file is part of Embed, a program implementing the Frozen Density Embedding (FDE) framework
! 
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at http://mozilla.org/MPL/2.0/.
!

module xml_parser
  use xml_structure
  implicit none

  private
  public xml_tag
  public xml_open,xml_close
  public xml_read,xml_write
  public xml_get_data,xml_add_data
  public xml_tag_open,xml_tag_close,xml_tag_next
  public xml_parent
  public set_attribute

  
  interface xml_tag_open
     module procedure tag_open
  end interface

  interface xml_tag_close
     module procedure tag_close
  end interface

  interface xml_tag_next
     module procedure next_tag,get_next_named
  end interface

  interface xml_add_data
     module procedure add_real_data, add_real_vector_data, add_char_data, add_char_data_vector
  end interface

  interface xml_get_data
     module procedure get_real_data, get_char_data  !,get_real_vector_data
  end interface

!radovan: this does not work with g95
! interface xml_attribute
!    module procedure get_attribute,set_attribute
! end interface

!radovan: this does not work with g95
! interface xml_name
!    module procedure get_name,set_name
! end interface
  
! A type for the enduser in which the data structure is hidden.
  
#include "xml_settings.h"

  type(xml_tag),pointer,save     :: xml_documents(:)
  type(xml_tag),pointer,save     :: xml_index
  integer,save               :: nodes_index,level


contains
  
  subroutine xml_init()
    nullify(xml_documents)
  end subroutine


  function tag_open(cur_tag,name,id)
    character(len=*) :: name
    character(len=*),optional :: id
    type(xml_tag),pointer :: cur_tag
    type(xml_tag),pointer :: tag_open

    if (.not.(associated(xml_documents))) cur_tag=xml_open('new')  
    tag_open=>add_child(cur_tag,name)
    if (present(id)) call add_attr(tag_open,'name',id)
  end function

  function tag_close(cur_tag)
    type(xml_tag),pointer :: cur_tag
    type(xml_tag),pointer :: tag_close
    if (associated(cur_tag)) then
      if(associated(cur_tag%parent)) then; tag_close=>cur_tag%parent    
      else; tag_close=>null(); end if
    end if
  end function
  
  subroutine set_attribute(cur_tag,name,value)
    type(xml_tag),pointer :: cur_tag
    character(len=*) :: name,value
    call add_attr(cur_tag,name,value)
  end subroutine

  character(len=attr_len) function get_attribute(cur_tag,name)
    type(xml_tag),pointer :: cur_tag
    character(len=*) :: name
    get_attribute=get_attr(cur_tag,name)
  end function

  subroutine set_name(cur_tag,name)
! Changes tagname
    type(xml_tag),pointer :: cur_tag
    character(len=*) :: name
    cur_tag%name=name
  end subroutine

  character(len=name_len) function get_name(cur_tag)
    type(xml_tag),pointer :: cur_tag
    get_name=cur_tag%name
  end function

  character(len=name_len) function xml_parent(cur_tag)
    type(xml_tag),pointer :: cur_tag
    if (associated(cur_tag%parent)) then
      xml_parent=cur_tag%parent%name
    else
      xml_parent=''
    end if
  end function

  subroutine add_real_data(cur_tag,real_data)
    type(xml_tag),pointer :: cur_tag
    real(kind=real_kind),intent(in) :: real_data(:,:)

    allocate(cur_tag%real_data(size(real_data,1),size(real_data,2)))
    cur_tag%real_data=real_data
    call add_attr(cur_tag,'type','real_data')
  end subroutine

  subroutine add_real_vector_data(cur_tag,real_data)
    type(xml_tag),pointer :: cur_tag
    real(kind=real_kind),intent(in) :: real_data(:)
    allocate(cur_tag%real_data(1,size(real_data)))
    cur_tag%real_data(1,:)=real_data
    call add_attr(cur_tag,'type','real_data')
  end subroutine

  subroutine add_char_data(cur_tag,char_data)
    type(xml_tag),pointer :: cur_tag
    character(len=*),intent(in) :: char_data
    call add_data(cur_tag,char_data)
  end subroutine

  subroutine add_char_data_vector(cur_tag,char_data)
    type(xml_tag),pointer :: cur_tag
    character(len=*),intent(in) :: char_data(:)
    integer :: i
    do i=1,size(char_data)
       call add_data(cur_tag,char_data(i))
    end do
  end subroutine


  subroutine get_real_data(tmp_data,cur_tag)
    type(xml_tag),intent(in) :: cur_tag
    real(kind=8),pointer :: tmp_data(:,:)
! We make a hard copy of the data into tmp_data.
! Free the memory tmp_data is pointing to.
    if (associated(tmp_data)) deallocate(tmp_data)
    tmp_data=>null()
    if (associated(cur_tag%real_data)) then 
       allocate(tmp_data(size(cur_tag%real_data,1),size(cur_tag%real_data,2)))
       tmp_data=cur_tag%real_data       
    end if
  end subroutine

!  function get_real_vector_data() result(tmp_data)
!    real(kind=real_kind),pointer :: tmp_data(:)
!
!    tmp_data=>null()
!    if (.not.(associated(nodes))) return
!    if (nodes_index>=size(nodes)) return
!    if (associated(nodes(nodes_index)%node%real_data)) tmp_data=>nodes(nodes_index)%node%real_data(1,:)       
!  end function

  subroutine get_char_data(tmp_data,cur_tag) 
    type(xml_tag),intent(in) :: cur_tag
    character(len=data_len),pointer :: tmp_data(:)
    tmp_data=>null()
    if (associated(cur_tag%data)) tmp_data=>cur_tag%data 
  end subroutine

  logical function next_tag(cur_tag,init)
    type(xml_tag),pointer :: cur_tag,tmp_tag
    integer,save     :: prev_index
    logical,optional :: init
    
    if (present(init)) level=0
    next_tag=.false.
    tmp_tag=>next_node(cur_tag,level)
    cur_tag=>tmp_tag
    if (associated(cur_tag)) then
      next_tag=.true.
    else
      if (level>0) level=0
      next_tag=.false.
    end if
  end function


  logical function get_next_named(cur_tag,name,attr,value)
    type(xml_tag),pointer :: cur_tag
    character(len=*) :: name
    character(len=*),optional :: attr,value
    logical :: add    

    get_next_named=.false.
    if (.not.associated(cur_tag)) return
    if(name=='') then;
       get_next_named=next_tag(cur_tag)
       return
    end if
    
    level=0
    do while(get_next_named.eqv..false.)
      add=.false.
      cur_tag=>next_node(cur_tag)
      if (.not.(associated(cur_tag))) return

      if (cur_tag%name==name) then
         if (present(attr).and.present(value)) then
            if (value==get_attr(cur_tag,attr)) add=.true.
         else; add=.true.
         end if
      end if

      if (add) then
        get_next_named=.true.
        exit
      end if
    end do
  end function



  function xml_close(doc_name)
     character(len=*),optional :: doc_name
     type(xml_tag),pointer         :: xml_close
     type(xml_tag),pointer         :: tmp_document
     type(xml_tag),pointer         :: tmp_xml_documents(:)
     integer :: i
     
     if (present(doc_name)) then
        do i=1,size(xml_documents)
           if (xml_documents(i)%name==doc_name) exit
        end do
        if (.not.(xml_documents(i)%name==doc_name)) then
           write(*,*) 'Warning: document name not a valid XML tree: ',doc_name
           return
        end if
     else; i=size(xml_documents) 
     end if
     tmp_document=>xml_documents(i)
     write(*,*) 'Deleting children'
     tmp_document=>del_children(tmp_document)
     write(*,*) 'Deleting documents'
     if (size(xml_documents)>1) then
	if (present(doc_name)) then
           do i=1,size(xml_documents)	
              if(doc_name==xml_documents(i)%name) exit
           end do
        else; i=size(xml_documents)
        end if
        allocate(tmp_xml_documents(size(xml_documents)-1))
        tmp_xml_documents(1:i-1)=xml_documents(1:i-1)
        tmp_xml_documents(i:)=xml_documents(i+1:)
        xml_documents=>tmp_xml_documents
     else
        deallocate(xml_documents)
        xml_documents=>null()
     end if
     xml_close=>null()
  end function

 
  function xml_open(doc_name)
    type(xml_tag),pointer :: xml_open
    character(len=*)  :: doc_name
    type(xml_tag),pointer :: tmp_document(:)
    integer           :: new_size

! Instantiate a new XML document tree.
    new_size=1
    if (associated(xml_documents)) then
       new_size=size(xml_documents)+1       
       allocate(tmp_document(new_size))
       tmp_document=xml_documents
       deallocate(xml_documents)
       xml_documents=>tmp_document
    else; allocate(xml_documents(new_size))
    end if        
! Name the document
    
    xml_documents(new_size)%name=doc_name
!    xml_open=>null()
    xml_open=>xml_documents(new_size)
    xml_open%child_index=0
    nullify(xml_open%data,xml_open%real_data,xml_open%attr,xml_open%children,xml_open%parent)
  end function

  function xml_read(file_name,doc_name)
    use xml_file

    character(len=*)         :: file_name
    character(len=*),optional:: doc_name
    character(len=name_len)  :: var
    character(len=attr_len)  :: attr
    character(len=1)         :: chr
    type(xml_tag),pointer    :: xml_read
    type(xml_tag),pointer    :: node
    integer                  :: ios=0,file_num,i,j,align
    integer(kind=8)          :: data_size, data_width,new_size
    logical                  :: is_tag,new_node
    type(xml_tag),pointer        :: tmp_xml_documents(:)
    
    is_tag=.false.
    new_node=.false.
    if(present(doc_name)) then; xml_read=>xml_open(doc_name)
    else; xml_read=>xml_open(file_name)
    end if

    node=>xml_read
! Point to the newly created document tree

    file_num=open_file(file_name)

    if (file_num<0) return
    do while(ios==0.or.ios==-2)
       var=''
       ios=read_var(file_num,var)
       select case(var)
       case ('<')
          is_tag=.true.
          ios=read_var(file_num,var)
          new_node=.true.
       case ('>')
          is_tag=.false.
          if (node%name=='dataset' .or. get_attr(node,'type')=='data') then
          ! Read a complet chunk of real numbers.
	     attr=get_attr(node,'size')
             read(attr, '(i8)') data_size
	     attr=get_attr(node,'width')
             read(attr, '(i8)') data_width
             if(data_width==0) data_width=1
             allocate(node%real_data(data_width,data_size))
! File pointer to next line
             read(file_num,'(a)',advance='yes') var

!            determine tab level
             chr=' '
             i=-1
             do while (chr==' ')
                read(file_num,'(1a)',advance='no') chr
                i=i+1
             end do
             backspace(file_num)

             if (i>0.and.(.not.(chr=='-'))) then; align=i-1
             else; align=i; end if             

             do i=1,data_size
                do j=1,align
                   read(file_num,'(a)',advance='no') chr
                end do
                do j=1,data_width-1
                   read(file_num,real_fmt,advance='no') node%real_data(j,i)
                end do
                read(file_num,real_fmt,advance='yes') node%real_data(data_width,i)
             end do
          end if
        case(' ')
           i=0
        case default
          if (is_tag) then          
             ios=read_var(file_num,attr)
             if (.not.(attr)=='=') then 
                write(*,*) 'ERROR: attribute not of the form: name="value"'
                write(*,*) 'TAG:      ', node%name
                write(*,*) 'ATTRIBUTE:', var
             end if
             ios=read_var(file_num,attr)
             call add_attr(node,var,attr)        
          else
             if (.not.(trim(var))=='') call add_data(node,trim(var))
          end if
       end select

       if (new_node) then
          if (var(1:1)=='/') then
            if (var(2:)==node%name(:)) then
               if (associated(node%parent)) then
                  node=>node%parent
                end if
            else
               write(*,*) "ERROR: Opening and closing tags are different"
               write(*,*) "Opening tag: ", node%name
               write(*,*) "Closing tag: ", var
               node=>xml_documents(size(xml_documents))
               node=>del_children(node)               
               return
            end if
          else
             node=>add_child(node,var)         
          end if
          new_node=.false.
       end if
    end do
    close(file_num)
  end function
  
  recursive subroutine write_child(file_num,node)
    use xml_file

    type(xml_tag),pointer :: node,tmp_node
    integer :: file_num
    integer,save :: tab_level=0
    character(len=attr_len) :: fmt

    integer :: i,j



    if(tab_level>0) then;
      write(fmt,'("(t",i4,",a)")') tab_level
    else; fmt='(a)'; end if;

    write(file_num,fmt,advance='no') '<'//trim(node%name)
    if(associated(node%attr)) then
       do i=1,size(node%attr,2)
          write(file_num, '(a)',advance='no') ' '//trim(node%attr(1,i))//'="'//trim(node%attr(2,i))//'"'
       end do
    end if
    write(file_num,'(a)',advance='yes') '>'

    tab_level=tab_level+tab_default
    if(tab_level>0) then;
      write(fmt,'("(t",i4,",a)")') tab_level
    else; fmt='(a)'; end if;

    if (associated(node%children)) then
       do i=1,size(node%children)
          tmp_node=>node%children(i)
          call write_child(file_num,tmp_node)
       end do
    end if

    if (associated(node%real_data)) then
       do i=1,size(node%real_data,2) 
          write(file_num,fmt,advance='no')
          do j=1,size(node%real_data,1)-1
             write(file_num,real_fmt,advance='no') node%real_data(j,i)
          end do
          write(file_num,real_fmt,advance='yes') node%real_data(j,i)
       end do
    end if
    if (associated(node%data)) then
       write(file_num,fmt,advance='no')
       do i=1,size(node%data)
          write(file_num,'(a)',advance='no') trim(node%data(i))
       end do
       write(file_num,fmt,advance='yes')
    end if


    tab_level=tab_level-tab_default
    if(tab_level>0) then;
      write(fmt,'("(t",i4,",a)")') tab_level
    else; fmt='(a)'; end if;

    write(file_num,fmt,advance='yes') '</'//trim(node%name)//'>'
    

  end subroutine


  subroutine xml_write(file_name,doc_name)
    use xml_file

    character(len=*)  :: file_name
    character(len=*),optional  :: doc_name
    integer :: file_num, doc_num,i
    type(xml_tag), pointer :: tmp_doc,tmp_node
    type(xml_tag),pointer :: tmp
    
  
    file_num=open_file(file_name)
    if (file_num<0) then;
       write(*,*) 'Failed to open file.'
       return
    end if
    if (present(doc_name)) then
       do doc_num=1,size(xml_documents)
          write(*,*) 'name:', xml_documents(doc_num)%name
          if (xml_documents(doc_num)%name==doc_name) exit
       end do
       if (.not.(xml_documents(doc_num)%name==doc_name)) then
          write(*,*) 'No document with name: ', doc_name
          return
       end if
    else; doc_num=size(xml_documents)
    end if 
    tmp_doc=>xml_documents(doc_num)
    do i=1,size(tmp_doc%children)
       tmp_node=>tmp_doc%children(i)
       call write_child(file_num,tmp_node)
    end do
    close(file_num)

  end subroutine
  


end module
