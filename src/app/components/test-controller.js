
class TestController {
  constructor() {
  	console.log("TestController")
  	var self = this;
  	this.message="HelloTedddst";
  	 $('#datetimepicker6').datetimepicker();
        $('#datetimepicker7').datetimepicker({
            useCurrent: false //Important! See issue #1075
        });
        $("#datetimepicker6").on("dp.change", function (e) {
            $('#datetimepicker7').data("DateTimePicker").minDate(e.date);
        });
        $("#datetimepicker7").on("dp.change", function (e) {
            $('#datetimepicker6').data("DateTimePicker").maxDate(e.date);
        });

   
  }
}
angular
  .module('app')
  .component('test', {
    templateUrl: 'app/containers/test.html',
    controller: TestController,
    controllerAs: "vm"
  })
 
