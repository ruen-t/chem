angular
  .module('app', ['ui.router', 'ngMaterial', 'angular-loading-bar', 'ui.sortable'])
  .controller('IndexController', IndexController)
  .directive("dropzone", function () {
    return {
      restrict: "A",
      link: function (scope, elem) {
        console.log(elem)
        elem.bind('dragover', function (e) {
          e.stopPropagation();
          e.preventDefault();
        });
        elem.bind('dragleave', function (e) {
          e.stopPropagation();
          e.preventDefault();
          scope.divClass = '';
        });
        elem.bind('drop', function (e) {
          e.stopPropagation();
          e.preventDefault();
          e.dataTransfer = e.originalEvent.dataTransfer;
          var files = e.dataTransfer.files;
          var file = files[0],
            reader = new FileReader();
          reader.onload = function (event) {
            var data = event.target.result;
            scope.vm.onDropFile(data)

            //holder.style.background = 'url(' + event.target.result + ') no-repeat center';
          };
          reader.readAsText(file);

        });
      }
    }
  })

function IndexController($state) {
  this.goto = $state.go

}
